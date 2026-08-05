#!/usr/bin/env python3
"""Append GAMBA-compatible phyloP baseline features to evaluation parquets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyBigWig

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import PHYLOP_NAMES, context_bounds, with_phylop

SUMMARY_NAMES = PHYLOP_NAMES


def read_scores(bigwig, chrom: str, start: int, end: int) -> np.ndarray:
    values = np.zeros(end - start, dtype=np.float64)
    try:
        intervals = bigwig.intervals(chrom, start, end)
    except RuntimeError:
        intervals = None
    if intervals is not None:
        for interval_start, interval_end, value in intervals:
            values[interval_start - start:interval_end - start] = value
    return np.asarray(np.round(values, 2), dtype=np.float32)


def read_oriented_scores(bigwig, row: dict, start: int,
                         end: int) -> np.ndarray:
    scores = read_scores(bigwig, row["chrom"], start, end)
    return scores[::-1] if row["strand"] == "-" else scores


def summarize(scores: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=np.float32)
    if scores.size == 0 or np.isnan(scores).all():
        return np.full(6, np.nan, dtype=np.float32)
    present = ~np.isnan(scores)
    positive = scores[(scores > 0) & present]
    negative = scores[(scores < 0) & present]
    denominator = float(present.sum())
    return np.asarray([
        np.nanmean(scores),
        np.nanstd(scores),
        float((scores > 0).sum()) / denominator if denominator else 0.0,
        float((scores < 0).sum()) / denominator if denominator else 0.0,
        np.nanmean(positive) if positive.size else 0.0,
        np.nanmean(negative) if negative.size else 0.0,
    ], dtype=np.float32)


def oriented_pool_span(row: dict) -> tuple[int, int]:
    if row["strand"] == "-":
        return (
            row["context_end"] - row["pool_end"],
            row["context_end"] - row["pool_start"],
        )
    return (
        row["context_start"] + row["pool_start"],
        row["context_start"] + row["pool_end"],
    )


def symmetric_context(row: dict, chrom_length: int,
                      length: int = 2048) -> tuple[int, int]:
    return context_bounds(
        chrom_length,
        row["start"],
        row["end"],
        row.get("strand", "+"),
        length,
        "symmetric",
    )


def baseline_roi_span(row: dict, chrom_length: int) -> tuple[int, int]:
    if row["scope"] == "100bp":
        return oriented_pool_span(row)
    if row["scope"] not in {"full", "roi"}:
        raise ValueError(f"unsupported scope: {row['scope']}")
    context_start, context_end = symmetric_context(row, chrom_length)
    start, end = sorted((row["start"], row["end"]))
    return max(start, context_start), min(end, context_end)


def annotate(input_path: Path, output_path: Path, bigwig_path: Path,
             batch_size: int) -> int:
    if not bigwig_path.exists():
        raise FileNotFoundError(
            f"{bigwig_path} is missing; run scripts/download_data.py first "
            "or pass --skip-phylop to the dataset processor"
        )
    source = pq.ParquetFile(input_path)
    for name in (
        *(f"phylop_{name}" for name in SUMMARY_NAMES),
        *(f"phylop_context_{name}" for name in SUMMARY_NAMES),
    ):
        if name in source.schema_arrow.names:
            raise ValueError(f"{input_path} already contains {name}")

    schema = with_phylop(source.schema_arrow)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    writer = pq.ParquetWriter(temporary, schema, compression="zstd")
    count = 0
    bigwig = pyBigWig.open(str(bigwig_path))
    chrom_lengths = bigwig.chroms()
    try:
        for batch in source.iter_batches(batch_size=batch_size):
            table = pa.Table.from_batches([batch])
            selected = table.select([
                "scope", "chrom", "start", "end", "strand", "context_start",
                "context_end", "pool_start", "pool_end",
            ]).to_pylist()
            pooled = np.empty((len(selected), 6), dtype=np.float32)
            context = np.empty((len(selected), 6), dtype=np.float32)
            for index, row in enumerate(selected):
                if row["chrom"] not in chrom_lengths:
                    raise ValueError(
                        f"{row['chrom']} is absent from {bigwig_path}"
                    )
                start, end = baseline_roi_span(
                    row, chrom_lengths[row["chrom"]]
                )
                pooled[index] = summarize(
                    read_oriented_scores(bigwig, row, start, end)
                )
                start, end = symmetric_context(row, chrom_lengths[row["chrom"]])
                context[index] = summarize(
                    read_oriented_scores(bigwig, row, start, end)
                )
            for index, name in enumerate(SUMMARY_NAMES):
                table = table.append_column(
                    f"phylop_{name}", pa.array(pooled[:, index], type=pa.float32())
                )
            for index, name in enumerate(SUMMARY_NAMES):
                table = table.append_column(
                    f"phylop_context_{name}",
                    pa.array(context[:, index], type=pa.float32()),
                )
            writer.write_table(table)
            count += len(selected)
    finally:
        bigwig.close()
        writer.close()
    temporary.replace(output_path)
    return count


def self_test() -> None:
    scores = np.asarray([0.0, 1.0, -2.0, 1.0], dtype=np.float32)
    expected = np.asarray([0.0, np.sqrt(1.5), 0.5, 0.25, 1.0, -2.0])
    np.testing.assert_allclose(summarize(scores), expected, rtol=1e-6)
    minus = {
        "strand": "-", "context_end": 1200, "pool_start": 100, "pool_end": 200
    }
    assert oriented_pool_span(minus) == (1000, 1100)
    assert symmetric_context({"start": 5000, "end": 5100}, 10_000) == (4026, 6074)
    long = {
        "scope": "full", "start": 4000, "end": 7000, "strand": "+",
    }
    assert baseline_roi_span(long, 10_000) == (4476, 6524)
    print("phyloP self-test passed")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument(
        "--bigwig",
        type=Path,
        default=root / "data/241-mammalian-2020v2.bigWig",
    )
    parser.add_argument("--output-dir", type=Path, default=root / "phylop")
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Replace each input parquet instead of writing under --output-dir.",
    )
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test()
    if not args.inputs:
        if args.self_test:
            return
        parser.error("provide at least one parquet input")
    if not args.bigwig.exists():
        raise FileNotFoundError(
            f"{args.bigwig} is missing; run scripts/download_data.py first"
        )
    for input_path in args.inputs:
        output = (
            input_path
            if args.in_place
            else args.output_dir / f"{input_path.stem}-phylop.parquet"
        )
        count = annotate(input_path, output, args.bigwig, args.batch_size)
        print(f"{output}: {count:,} rows")


if __name__ == "__main__":
    main()
