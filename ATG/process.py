#!/usr/bin/env python3
"""Minimal GAMBA ATG-5-way source and parquet builder.

Label construction and chromosome-balanced sampling are derived from Microsoft
GAMBA's make_ATG_data.py and ATG_reps.py at commit e83984e.
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ATG_SCHEMA, COMPLEMENT, context_for_region, write_parquet
from scripts.add_phylop import annotate


LABELS = {
    1: ("label1_start_pos", None, "start"),
    2: ("label2_noncoding_near_pos", "label2_delta_bp", "noncoding_near"),
    3: ("label3_noncoding_far_pos", "label3_delta_bp", "noncoding_far"),
    4: ("label4_same_inframe_met_pos", "label4_delta_bp", "inframe_methionine"),
    5: ("label5_same_outframe_atg_pos", "label5_delta_bp", "outframe_atg"),
}
CHROMS = [f"chr{i}" for i in range(1, 23)]


def attrs(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        item = item.strip()
        if " " in item:
            key, value = item.split(" ", 1)
            result[key] = value.strip().strip('"')
    return result


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    strand: str
    cds: list[tuple[int, int, int]] = field(default_factory=list)
    blocks: list[tuple[int, int, int]] = field(default_factory=list)

    def finish(self) -> None:
        segments = sorted(self.cds, reverse=self.strand == "-")
        phase = segments[0][2] if segments and segments[0][2] in (0, 1, 2) else 0
        offset = 0
        for index, (start, end, _) in enumerate(segments):
            if index == 0 and phase:
                if self.strand == "+":
                    start = min(end, start + phase)
                else:
                    end = max(start, end - phase)
            if start < end:
                self.blocks.append((start, end, offset))
                offset += end - start

    def sequence(self, chromosome: str) -> str:
        sequence = "".join(chromosome[start:end] for start, end, _ in self.blocks)
        return sequence.translate(COMPLEMENT)[::-1] if self.strand == "-" else sequence

    def genomic_position(self, offset: int) -> int | None:
        for start, end, block_offset in self.blocks:
            if block_offset <= offset < block_offset + end - start:
                distance = offset - block_offset
                return start + distance if self.strand == "+" else end - 1 - distance
        return None


def transcripts(path: Path, chrom: str) -> list[Transcript]:
    result = {}
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if (
                len(fields) < 9 or fields[0] != chrom or fields[2] != "CDS"
                or fields[6] not in {"+", "-"} or 'tag "MANE_Select"' not in fields[8]
            ):
                continue
            attributes = attrs(fields[8])
            transcript_id, gene_id = attributes.get("transcript_id"), attributes.get("gene_id")
            if not transcript_id or not gene_id:
                continue
            transcript = result.setdefault(
                transcript_id,
                Transcript(transcript_id, gene_id, fields[6]),
            )
            try:
                phase = int(fields[7])
            except ValueError:
                phase = 0
            transcript.cds.append((int(fields[3]) - 1, int(fields[4]), phase))
    for transcript in result.values():
        transcript.finish()
    return [transcript for transcript in result.values()
            if sum(end - start for start, end, _ in transcript.blocks) >= 3]


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def outside_cds_atgs(chromosome: str, intervals: list[tuple[int, int]]) -> list[int]:
    starts = [start for start, _ in intervals]
    positions, index = [], chromosome.find("ATG")
    while index != -1:
        interval = bisect_right(starts, index) - 1
        if interval < 0 or index + 3 > intervals[interval][1]:
            positions.append(index)
        index = chromosome.find("ATG", index + 1)
    return positions


def closest_at_distance(position: int, values: list[int], minimum: int,
                        maximum: int | None = None) -> tuple[int, int] | None:
    candidates = []
    right = bisect_left(values, position + minimum)
    if right < len(values):
        delta = values[right] - position
        if maximum is None or delta <= maximum:
            candidates.append((values[right], delta))
    left = bisect_right(values, position - minimum) - 1
    if left >= 0:
        delta = position - values[left]
        if maximum is None or delta <= maximum:
            candidates.append((values[left], -delta))
    return min(candidates, key=lambda item: abs(item[1])) if candidates else None


def transcript_sites(transcript: Transcript, chromosome: str) -> tuple[int, list[int], list[int]] | None:
    sequence = transcript.sequence(chromosome).upper()
    if not sequence.startswith("ATG"):
        return None
    start = transcript.genomic_position(0)
    if start is None:
        return None
    inframe, outframe = [], []
    for offset in range(len(sequence) - 2):
        if sequence[offset:offset + 3] != "ATG":
            continue
        position = transcript.genomic_position(offset)
        if position is None:
            continue
        if offset % 3 == 0:
            if offset:
                inframe.append(position)
        else:
            outframe.append(position)
    return start, sorted(inframe), sorted(outframe)


def generate_chromosome(chrom: str, gtf: Path, genome: Fasta, output: Path,
                        sample_size: int, seed: int) -> int:
    chromosome = genome[chrom][:].seq.upper()
    chrom_transcripts = transcripts(gtf, chrom)
    cds = merge([
        (start, end)
        for transcript in chrom_transcripts
        for start, end, _ in transcript.blocks
    ])
    noncoding = outside_cds_atgs(chromosome, cds)
    sites = []
    for transcript in chrom_transcripts:
        found = transcript_sites(transcript, chromosome)
        if found:
            sites.append((transcript, *found))
    if sample_size and sample_size < len(sites):
        sites = random.Random(seed).sample(sites, sample_size)

    output.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "chrom", "transcript_id", "gene_id", "strand", "label1_start_pos",
        "label2_noncoding_near_pos", "label2_delta_bp",
        "label3_noncoding_far_pos", "label3_delta_bp",
        "label4_same_inframe_met_pos", "label4_delta_bp",
        "label5_same_outframe_atg_pos", "label5_delta_bp",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for transcript, anchor, inframe, outframe in sites:
            near = closest_at_distance(anchor, noncoding, 2000, 5000)
            far = closest_at_distance(anchor, noncoding, 100_000)
            same_frame = (
                min(((position, position - anchor) for position in inframe),
                    key=lambda item: abs(item[1]))
                if inframe else None
            )
            other_frame = (
                min(((position, position - anchor) for position in outframe),
                    key=lambda item: abs(item[1]))
                if outframe else None
            )
            values = []
            for item in (near, far, same_frame, other_frame):
                values.extend(item if item else (".", "."))
            writer.writerow([
                chrom, transcript.transcript_id, transcript.gene_id,
                transcript.strand, anchor, *values,
            ])
    return len(sites)


def concatenate(paths: list[Path], output: Path) -> None:
    with output.open("w", newline="") as target:
        for index, path in enumerate(sorted(paths)):
            with path.open(newline="") as source:
                if index:
                    next(source)
                for line in source:
                    target.write(line)


def even_sample(frame: pd.DataFrame, total: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    frame = frame[frame.chrom.isin(CHROMS)].copy()
    base, remainder = divmod(total, len(CHROMS))
    targets = {chrom: base + (index < remainder) for index, chrom in enumerate(CHROMS)}
    selected, remaining = [], frame.copy()
    while any(targets.values()):
        carry = 0
        for chrom in CHROMS:
            wanted = targets[chrom]
            if not wanted:
                continue
            available = remaining[remaining.chrom == chrom]
            take = min(wanted, len(available))
            if take:
                indices = rng.choice(available.index.to_numpy(), size=take, replace=False)
                selected.append(remaining.loc[indices])
                remaining = remaining.drop(indices)
            carry += wanted - take
            targets[chrom] = 0
        if not carry or remaining.empty:
            break
        available_chroms = [chrom for chrom in CHROMS if (remaining.chrom == chrom).any()]
        for index in range(carry):
            targets[available_chroms[index % len(available_chroms)]] += 1
    result = pd.concat(selected) if selected else frame.iloc[:0]
    return result.sample(total, random_state=seed) if len(result) > total else result


def parquet_rows(
    sampled: pd.DataFrame,
    genome: Fasta,
    context_policy: str,
):
    for row in sampled.itertuples(index=False):
        positions = {label_id: int(getattr(row, column))
                     for label_id, (column, _, _) in LABELS.items()}
        pair_id = f"{row.chrom}|{row.transcript_id}|{row.strand}|{positions[1]}"
        contexts = []
        for label_id, position in positions.items():
            context = context_for_region(
                genome,
                row.chrom,
                position,
                position + 3,
                row.strand,
                policy=context_policy,
            )
            if context is None:
                contexts = []
                break
            contexts.append((label_id, position, context))
        for label_id, position, context in contexts:
            _, delta_column, label = LABELS[label_id]
            yield {
                "split": "test",
                "sequence": context["sequence"],
                "label": label,
                "pair_id": pair_id,
                "category": "ATG",
                "scope": "roi",
                "context_policy": context_policy,
                "chrom": row.chrom,
                "start": position,
                "end": position + 3,
                "strand": row.strand,
                "context_start": context["context_start"],
                "context_end": context["context_end"],
                "roi_start": context["roi_start"],
                "roi_end": context["roi_end"],
                "pool_start": context["roi_start"],
                "pool_end": context["roi_end"],
                "name": f"{row.transcript_id}_L{label_id}",
                "transcript_id": row.transcript_id,
                "gene_id": row.gene_id,
                "label_id": label_id,
                "delta_bp": 0 if delta_column is None else int(getattr(row, delta_column)),
            }


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=root.parent / "data")
    parser.add_argument("--source-dir", type=Path, default=root / "source")
    parser.add_argument("--output-dir", type=Path, default=root)
    parser.add_argument("--n-examples", type=int, default=2000)
    parser.add_argument("--per-chrom-limit", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-source-generation", action="store_true")
    parser.add_argument(
        "--bigwig",
        type=Path,
        default=root.parent / "data/241-mammalian-2020v2.bigWig",
    )
    parser.add_argument("--skip-phylop", action="store_true")
    args = parser.parse_args()

    genome = Fasta(str(args.data_root / "hg38.ml.fa"))
    args.source_dir.mkdir(parents=True, exist_ok=True)
    chromosome_paths = [
        args.source_dir / f"{chrom}_atg_5way_labels.tsv" for chrom in CHROMS
    ]
    if not args.skip_source_generation:
        for chrom, path in zip(CHROMS, chromosome_paths):
            count = generate_chromosome(
                chrom, args.data_root / "gtfs" / f"{chrom}.gtf",
                genome, path, args.per_chrom_limit, args.seed,
            )
            print(f"{chrom}: {count:,} source rows")

    combined = args.source_dir / "all_chr_atg_5way.tsv"
    concatenate(chromosome_paths, combined)
    frame = pd.read_csv(combined, sep="\t")
    for column, _, _ in LABELS.values():
        frame = frame[frame[column].astype(str) != "."]
    sampled = even_sample(frame, args.n_examples, args.seed)
    sampled.to_csv(args.source_dir / "sampled_examples_atg5.tsv", sep="\t", index=False)
    for context_policy, filename_policy in (
        ("causal", "causal"),
        ("symmetric", "bidi"),
    ):
        output = args.output_dir / f"atg-gamba-{filename_policy}.parquet"
        count = write_parquet(
            output,
            parquet_rows(sampled, genome, context_policy),
            ATG_SCHEMA,
        )
        if not args.skip_phylop:
            annotate(output, output, args.bigwig, 2048)
        assert count == 5 * len(sampled)
        print(
            f"ATG-5-way ({context_policy}): "
            f"{len(sampled):,} examples, {count:,} rows"
        )
    genome.close()


if __name__ == "__main__":
    main()
