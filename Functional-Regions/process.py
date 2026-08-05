#!/usr/bin/env python3
"""Minimal GAMBA functional-region generator and parquet builder.

Region construction is derived from Microsoft GAMBA's
data_processing/create_eval_data.py at commit e83984e.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import shutil
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import context_for_region, stable_100bp_span, write_parquet
from scripts.add_phylop import annotate


GAMBA_CATEGORIES = [
    "repeats", "UCNE", "vista_enhancer", "promoters", "UTR5", "UTR3",
    "coding_regions", "exons", "introns", "upstream_TSS",
]
TEST_CHROMS = {"chr2", "chr3", "chr16", "chr22"}
ALL_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
ROLES = {
    "feature": "{category}",
    "upstream": "{category}_upstream",
    "random": "{category}_random",
    "random-noannot": "{category}_random-noannot",
}


def overlaps(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    index = bisect_left(intervals, (start, end))
    return (
        (index > 0 and intervals[index - 1][1] > start)
        or (index < len(intervals) and intervals[index][0] < end)
    )


def insert_merge(intervals: list[tuple[int, int]], start: int, end: int) -> None:
    index = bisect_left(intervals, (start, end))
    left = index - 1
    while left >= 0 and intervals[left][1] > start:
        start, end = min(start, intervals[left][0]), max(end, intervals[left][1])
        left -= 1
    left += 1
    right = index
    while right < len(intervals) and intervals[right][0] < end:
        start, end = min(start, intervals[right][0]), max(end, intervals[right][1])
        right += 1
    intervals[left:right] = [(start, end)]


def normalize_chrom(chrom: str, canonical: set[str]) -> str:
    if chrom in canonical:
        return chrom
    alternate = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if alternate in canonical:
        return alternate
    if chrom in {"M", "MT", "chrM", "chrMT"}:
        return next((name for name in ("chrM", "MT", "M") if name in canonical), chrom)
    return chrom


def bed_rows(path: Path, category: str, canonical: set[str],
             chromosomes: set[str], blocked: set[str] = set()) -> list[dict]:
    rows = []
    with path.open() as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields or fields[0].startswith("#") or len(fields) < 3:
                continue
            chrom = normalize_chrom(fields[0], canonical)
            name = fields[3] if len(fields) > 3 else category
            if chrom not in chromosomes or name in blocked:
                continue
            start, end = int(fields[1]), int(fields[2])
            if end < start:
                start, end = end, start
            if end == start:
                continue
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "score": float(fields[4]) if len(fields) > 4 else 0.0,
                "strand": fields[5] if len(fields) > 5 else ".",
                "category": category,
            })
    return rows


def blocked_ucnes(path: Path) -> set[str]:
    blocked = set()
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            group = line.replace(",", " ").split()
            blocked.update(group[1:])
    return blocked


def vista_rows(path: Path, canonical: set[str], chromosomes: set[str]) -> list[dict]:
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        coordinate = next(
            (name for name in ("coordinate_hg38", "coord", "Element Coordinates")
             if name in reader.fieldnames),
            None,
        )
        if coordinate is None:
            raise ValueError(f"{path}: missing a VISTA coordinate column")
        identifier = "vista_id" if "vista_id" in reader.fieldnames else "Element ID"
        for row in reader:
            value = row.get(coordinate)
            if not value:
                continue
            chrom, positions = value.strip().split(":")
            start, end = map(int, positions.split("-"))
            chrom = normalize_chrom(chrom, canonical)
            if chrom not in chromosomes:
                continue
            rows.append({
                "chrom": chrom,
                "start": start - 1,
                "end": end,
                "name": row[identifier],
                "score": 0.0,
                "strand": row.get("strand") or ".",
                "category": "vista_enhancer",
            })
    return rows


def attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        item = item.strip()
        if " " in item:
            key, value = item.split(" ", 1)
            result[key] = value.strip().strip('"')
    return result


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def gtf_features(gtf_dir: Path, chromosomes: set[str]) -> tuple[dict[str, list[dict]], dict[str, list[tuple[int, int]]]]:
    features = {name: [] for name in GAMBA_CATEGORIES[-4:]}
    annotated = defaultdict(list)
    paths = sorted(gtf_dir.glob("*.gtf")) + sorted(gtf_dir.glob("*.gtf.gz"))
    for path in paths:
        if path.name.split(".")[0] not in chromosomes:
            continue
        transcripts = {}
        with open_text(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, kind, start, end, _, strand, _, attrs_text = fields
                attrs = attributes(attrs_text)
                transcript_id = attrs.get("transcript_id")
                if kind == "transcript" and transcript_id and attrs.get("gene_id") and attrs.get("transcript_type"):
                    transcripts[transcript_id] = {
                        "chrom": chrom,
                        "start": int(start) - 1,
                        "end": int(end),
                        "strand": strand,
                        "type": attrs["transcript_type"],
                        "canonical": (
                            'tag "Ensembl_canonical"' in attrs_text
                            or 'tag "CCDS"' in attrs_text
                        ),
                        "exons": [],
                        "cds": [],
                    }
                elif kind in {"exon", "CDS"} and transcript_id in transcripts:
                    transcripts[transcript_id]["exons" if kind == "exon" else "cds"].append(
                        (int(start) - 1, int(end))
                    )

        for transcript_id, transcript in transcripts.items():
            if not transcript["canonical"]:
                continue
            chrom, strand = transcript["chrom"], transcript["strand"]
            base = {"chrom": chrom, "strand": strand, "name": transcript_id, "score": 0.0}
            if transcript["type"] == "protein_coding":
                for start, end in transcript["cds"]:
                    features["coding_regions"].append({**base, "start": start, "end": end})
            for start, end in transcript["exons"]:
                features["exons"].append({**base, "start": start, "end": end})
            exons = sorted(transcript["exons"])
            for (_, intron_start), (intron_end, _) in zip(exons, exons[1:]):
                if intron_end > intron_start:
                    features["introns"].append({
                        **base, "start": intron_start, "end": intron_end,
                    })
            if strand == "+":
                start, end = transcript["start"] - 2000, transcript["start"]
            else:
                start, end = transcript["end"], transcript["end"] + 2000
            if end > start:
                features["upstream_TSS"].append({**base, "start": start, "end": end})
            annotated[chrom].extend(transcript["cds"] + transcript["exons"])

    for category, rows in features.items():
        for row in rows:
            row["category"] = category
    return features, annotated


def noncoding_rows(annotated: dict[str, list[tuple[int, int]]],
                   chrom_lengths: dict[str, int]) -> list[dict]:
    rows = []
    for chrom, intervals in annotated.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        cursor = 0
        for start, end in merged:
            if start > cursor:
                rows.append({
                    "chrom": chrom, "start": cursor, "end": start,
                    "name": "noncoding_regions", "score": 0.0, "strand": ".",
                    "category": "noncoding_regions",
                })
            cursor = max(cursor, end)
        if cursor < chrom_lengths[chrom]:
            rows.append({
                "chrom": chrom, "start": cursor, "end": chrom_lengths[chrom],
                "name": "noncoding_regions", "score": 0.0, "strand": ".",
                "category": "noncoding_regions",
            })
    return rows


def retain_nonoverlapping(categories: dict[str, list[dict]], order: list[str],
                          limit: int, seed: int) -> dict[str, list[dict]]:
    rng, occupied, output = random.Random(seed), defaultdict(list), {}
    for category in order:
        pool = list(categories[category])
        rng.shuffle(pool)
        kept = []
        for row in pool:
            start, end = int(row["start"]), int(row["end"])
            if end <= start or overlaps(occupied[row["chrom"]], start, end):
                continue
            kept.append(row)
            insert_merge(occupied[row["chrom"]], start, end)
            if limit and len(kept) >= limit:
                break
        output[category] = kept
    return output


def add_upstream(categories: dict[str, list[dict]], chrom_lengths: dict[str, int],
                 distance: int) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    occupied = defaultdict(lambda: defaultdict(list))
    for category, rows in categories.items():
        for row in rows:
            insert_merge(occupied[category][row["chrom"]], row["start"], row["end"])

    anchors = {category: [] for category in categories}
    upstream = {f"{category}_upstream": [] for category in categories}
    pair_id = 0
    for category, rows in categories.items():
        for row in rows:
            length = row["end"] - row["start"]
            if row["strand"] == "-":
                start = row["end"] + distance
                end = start + length
            else:
                end = row["start"] - distance
                start = end - length
            chrom = row["chrom"]
            if (
                start < 0 or end > chrom_lengths[chrom]
                or overlaps(occupied[category][chrom], start, end)
            ):
                continue
            anchor = {**row, "pair_id": pair_id}
            anchors[category].append(anchor)
            upstream[f"{category}_upstream"].append({
                **row,
                "start": start,
                "end": end,
                "name": f"{row['name']}_up",
                "score": 0.0,
                "category": f"{category}_upstream",
                "pair_id": pair_id,
            })
            pair_id += 1
    return anchors, upstream


def random_controls(categories: dict[str, list[dict]], chrom_lengths: dict[str, int],
                    seed: int, attempts: int) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    all_occupied, category_occupied = defaultdict(list), defaultdict(lambda: defaultdict(list))
    for category, rows in categories.items():
        for row in rows:
            insert_merge(all_occupied[row["chrom"]], row["start"], row["end"])
            insert_merge(category_occupied[category][row["chrom"]], row["start"], row["end"])

    rng = random.Random(seed)
    noannot_occupied = defaultdict(list, {
        chrom: list(intervals) for chrom, intervals in all_occupied.items()
    })
    used_random = defaultdict(list)
    noannot = {f"{category}_random-noannot": [] for category in categories}
    randoms = {f"{category}_random": [] for category in categories}

    def sample(chrom: str, length: int, occupied: list[tuple[int, int]]) -> tuple[int, int] | None:
        for _ in range(attempts):
            start = rng.randint(0, chrom_lengths[chrom] - length)
            end = start + length
            if not overlaps(occupied, start, end):
                return start, end
        return None

    for category, rows in categories.items():
        category_random_occupied = defaultdict(list)
        for chrom, intervals in category_occupied[category].items():
            category_random_occupied[chrom] = list(intervals)
            for interval in used_random[chrom]:
                insert_merge(category_random_occupied[chrom], *interval)
        for row in rows:
            chrom, length = row["chrom"], row["end"] - row["start"]
            sampled = sample(chrom, length, noannot_occupied[chrom])
            if sampled:
                insert_merge(noannot_occupied[chrom], *sampled)
                noannot[f"{category}_random-noannot"].append({
                    **row, "start": sampled[0], "end": sampled[1],
                    "name": f"{row['name']}_rno", "score": 0.0,
                    "category": f"{category}_random-noannot",
                })

            combined = category_random_occupied[chrom]
            sampled = sample(chrom, length, combined)
            if sampled:
                insert_merge(combined, *sampled)
                insert_merge(used_random[chrom], *sampled)
                randoms[f"{category}_random"].append({
                    **row, "start": sampled[0], "end": sampled[1],
                    "name": f"{row['name']}_rcat", "score": 0.0,
                    "category": f"{category}_random",
                })
    return noannot, randoms


def write_beds(root: Path, groups: dict[str, list[dict]]) -> None:
    for category, rows in groups.items():
        by_chrom = defaultdict(list)
        for row in rows:
            by_chrom[row["chrom"]].append(row)
        for chrom, chrom_rows in by_chrom.items():
            path = root / category / f"{chrom}.bed"
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w") as handle:
                for row in chrom_rows:
                    handle.write("\t".join(map(str, [
                        row["chrom"], row["start"], row["end"], row["name"],
                        row["score"], row["strand"], row["pair_id"],
                    ])) + "\n")


def build_regions(data_root: Path, output: Path, include_noncoding: bool,
                  limit: int, seed: int) -> list[str]:
    genome = Fasta(str(data_root / "hg38.ml.fa"))
    canonical = set(genome.keys())
    chromosomes = set(canonical)
    chrom_lengths = {chrom: len(genome[chrom]) for chrom in canonical}
    region_info = data_root / "region_info"
    gtf, annotated = gtf_features(data_root / "gtfs", chromosomes)

    categories = {
        "repeats": bed_rows(region_info / "repeats_hg38.bed", "repeats", canonical, chromosomes),
        "UCNE": bed_rows(
            region_info / "hg38_UCNE_coordinates.bed", "UCNE", canonical, chromosomes,
            blocked_ucnes(region_info / "ucne_paralogues.txt"),
        ),
        "vista_enhancer": vista_rows(region_info / "experiments.tsv", canonical, chromosomes),
        "promoters": bed_rows(region_info / "promoters.bed", "promoters", canonical, chromosomes),
        "UTR5": bed_rows(region_info / "UCSC_5UTR_exons.bed", "UTR5", canonical, chromosomes),
        "UTR3": bed_rows(region_info / "UCSC_3UTR_exons.bed", "UTR3", canonical, chromosomes),
        **gtf,
    }
    order = list(GAMBA_CATEGORIES)
    if include_noncoding:
        categories["noncoding_regions"] = noncoding_rows(annotated, chrom_lengths)
        order.append("noncoding_regions")

    retained = retain_nonoverlapping(categories, order, limit, seed)
    if include_noncoding:
        canonical = {name: retained[name] for name in GAMBA_CATEGORIES}
        canonical, upstream = add_upstream(canonical, chrom_lengths, 2000)
        noannot, randoms = random_controls(canonical, chrom_lengths, seed, 2000)

        all_anchors, all_upstream = add_upstream(retained, chrom_lengths, 2000)
        all_noannot, all_randoms = random_controls(
            all_anchors, chrom_lengths, seed, 2000
        )
        retained = {**canonical, "noncoding_regions": all_anchors["noncoding_regions"]}
        upstream["noncoding_regions_upstream"] = all_upstream[
            "noncoding_regions_upstream"
        ]
        noannot["noncoding_regions_random-noannot"] = all_noannot[
            "noncoding_regions_random-noannot"
        ]
        randoms["noncoding_regions_random"] = all_randoms[
            "noncoding_regions_random"
        ]
    else:
        retained, upstream = add_upstream(retained, chrom_lengths, 2000)
        noannot, randoms = random_controls(retained, chrom_lengths, seed, 2000)
    shutil.rmtree(output, ignore_errors=True)
    output.mkdir(parents=True)
    write_beds(output, {**retained, **upstream, **noannot, **randoms})
    genome.close()
    return order


def read_role_maps(root: Path, category: str, chromosomes: set[str]) -> dict[str, dict[str, dict]]:
    maps = {role: {} for role in ROLES}
    for role, template in ROLES.items():
        directory = root / template.format(category=category)
        for path in sorted(directory.glob("chr*.bed")):
            if path.stem not in chromosomes:
                continue
            with path.open() as handle:
                for fields in csv.reader(handle, delimiter="\t"):
                    maps[role][fields[6]] = {
                        "chrom": fields[0], "start": int(fields[1]), "end": int(fields[2]),
                        "name": fields[3], "strand": fields[5], "pair_id": fields[6],
                    }
    return maps


def parquet_row(
    genome,
    category: str,
    label: str,
    scope: str,
    row: dict,
    context_policy: str,
    pool: tuple[int, int] | None = None,
) -> dict | None:
    context = context_for_region(
        genome,
        row["chrom"],
        row["start"],
        row["end"],
        row["strand"],
        policy=context_policy,
    )
    if context is None:
        return None
    pool = pool or (context["roi_start"], context["roi_end"])
    return {
        "split": "test" if row["chrom"] in TEST_CHROMS else "train",
        "sequence": context["sequence"],
        "label": label,
        "pair_id": row["pair_id"],
        "category": category,
        "scope": scope,
        "context_policy": context_policy,
        "chrom": row["chrom"],
        "start": row["start"],
        "end": row["end"],
        "strand": row["strand"],
        "context_start": context["context_start"],
        "context_end": context["context_end"],
        "roi_start": context["roi_start"],
        "roi_end": context["roi_end"],
        "pool_start": pool[0],
        "pool_end": pool[1],
        "name": row["name"],
    }


def build_parquets(regions: Path, genome_path: Path, output: Path,
                   categories: list[str], seed: int,
                   bigwig_path: Path | None) -> None:
    genome = Fasta(str(genome_path))
    chromosomes = set(ALL_CHROMS)
    def binary_rows(control: str, context_policy: str):
        for category in categories:
            maps = read_role_maps(regions, category, chromosomes)
            common = set.intersection(*(set(values) for values in maps.values()))
            for pair_id in sorted(common):
                for role, label in (("feature", "feature"), (control, control)):
                    result = parquet_row(
                        genome,
                        category,
                        label,
                        "full",
                        maps[role][pair_id],
                        context_policy,
                    )
                    if result:
                        yield result

    def multiclass_rows(scope: str, context_policy: str):
        for category in categories:
            maps = read_role_maps(regions, category, chromosomes)
            common = set.intersection(*(set(values) for values in maps.values()))
            for pair_id in sorted(common):
                row = maps["feature"][pair_id]
                result = parquet_row(
                    genome, category, category, scope, row, context_policy
                )
                if result and scope == "full":
                    yield result
                elif result:
                    span = stable_100bp_span(
                        category, pair_id, result["roi_start"], result["roi_end"], seed
                    )
                    if span:
                        yield {**result, "pool_start": span[0], "pool_end": span[1]}

    output.mkdir(parents=True, exist_ok=True)
    for context_policy, filename_policy in (
        ("causal", "causal"),
        ("symmetric", "bidi"),
    ):
        for control in ("upstream", "random", "random-noannot"):
            name = (
                f"functional-{control}-gamba-{filename_policy}.parquet"
            )
            count = write_parquet(
                output / name,
                binary_rows(control, context_policy),
            )
            if bigwig_path is not None:
                annotate(output / name, output / name, bigwig_path, 2048)
            print(f"{name}: {count:,} rows")
        for scope, label in (("full", "full"), ("100bp", "100bp")):
            name = (
                f"functional-multiclass-gamba-{label}-"
                f"{filename_policy}.parquet"
            )
            count = write_parquet(
                output / name, multiclass_rows(scope, context_policy)
            )
            if bigwig_path is not None:
                annotate(output / name, output / name, bigwig_path, 2048)
            print(f"{name}: {count:,} rows")
    genome.close()


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=root.parent / "data")
    parser.add_argument("--regions-dir", type=Path, default=root / "regions")
    parser.add_argument("--output-dir", type=Path, default=root)
    parser.add_argument("--limit-per-category", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--bigwig",
        type=Path,
        default=root.parent / "data/241-mammalian-2020v2.bigWig",
    )
    parser.add_argument("--skip-phylop", action="store_true")
    parser.add_argument(
        "--include-noncoding", action="store_true"
    )
    parser.add_argument("--skip-regions", action="store_true")
    parser.add_argument("--skip-parquets", action="store_true")
    args = parser.parse_args()

    categories = list(GAMBA_CATEGORIES)
    if not args.skip_regions:
        categories = build_regions(
            args.data_root, args.regions_dir, args.include_noncoding,
            args.limit_per_category, args.seed,
        )
    elif args.include_noncoding:
        categories.append("noncoding_regions")
    if not args.skip_parquets:
        build_parquets(
            args.regions_dir, args.data_root / "hg38.ml.fa", args.output_dir,
            categories, args.seed,
            None if args.skip_phylop else args.bigwig,
        )


if __name__ == "__main__":
    main()
