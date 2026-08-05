#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (
    ATG_SCHEMA,
    PHYLOP_FIELDS,
    SCHEMA,
    context_bounds,
    with_phylop,
)


def compare_trees(expected: Path, actual: Path) -> None:
    expected_files = sorted(path.relative_to(expected) for path in expected.rglob("*.bed"))
    actual_files = sorted(path.relative_to(actual) for path in actual.rglob("*.bed"))
    assert actual_files == expected_files
    mismatches = [
        path for path in expected_files
        if (expected / path).read_bytes() != (actual / path).read_bytes()
    ]
    assert not mismatches, f"BED mismatches: {mismatches[:10]}"


def verify_parquet(path: Path, require_phylop: bool) -> int:
    table = pq.read_table(path)
    expected = ATG_SCHEMA if path.name.startswith("atg-gamba-") else SCHEMA
    if require_phylop:
        expected = with_phylop(expected)
    assert table.schema == expected
    frame = table.to_pandas()
    filename_policy = path.stem.rsplit("-", 1)[-1]
    assert filename_policy in {"causal", "bidi"}
    context_policy = "symmetric" if filename_policy == "bidi" else "causal"
    assert set(frame.context_policy) == {context_policy}
    assert frame.sequence.str.len().between(1, 2048).all()
    assert (frame.roi_start >= 0).all()
    assert (frame.roi_end > frame.roi_start).all()
    assert (frame.pool_start >= 0).all()
    assert (frame.pool_end > frame.pool_start).all()
    assert (frame.pool_end <= frame.sequence.str.len()).all()
    assert (frame.pool_start >= frame.roi_start).all()
    assert (frame.pool_end <= frame.roi_end).all()
    assert not frame.duplicated(
        ["scope", "category", "label", "pair_id"]
    ).any()
    if require_phylop:
        assert not frame[[field.name for field in PHYLOP_FIELDS]].isna().any().any()

    if path.name.startswith("functional-") and "multiclass" not in path.name:
        assert set(frame.split) == {"train", "test"}
        assert set(frame.loc[frame.split == "test", "chrom"]) == {
            "chr2", "chr3", "chr16", "chr22"
        }
        assert "noncoding_regions" in set(frame.category)
        grouped = frame.groupby(["category", "pair_id"])
        groups = grouped.label.agg(set)
        control = path.stem.removeprefix("functional-").split("-gamba", 1)[0]
        expected = {"feature", control}
        assert groups.map(lambda labels: labels == expected).all()
        assert grouped.size().eq(2).all()
        assert len(frame[frame.category != "noncoding_regions"]) == 169_674
    elif "multiclass" in path.name:
        assert set(frame.split) == {"train", "test"}
        assert set(frame.loc[frame.split == "test", "chrom"]) == {
            "chr2", "chr3", "chr16", "chr22"
        }
        assert "noncoding_regions" in set(frame.category)
        expected_scope = "100bp" if "-100bp" in path.name else "full"
        assert set(frame.scope) == {expected_scope}
        assert frame.label.eq(frame.category).all()
        if expected_scope == "100bp":
            assert (frame.pool_end - frame.pool_start).eq(100).all()
            assert len(frame[frame.category != "noncoding_regions"]) == 58_643
        else:
            assert len(frame[frame.category != "noncoding_regions"]) == 84_837
    else:
        grouped = frame.groupby("pair_id")
        groups = grouped.label.agg(set)
        expected = {
            "start", "noncoding_near", "noncoding_far",
            "inframe_methionine", "outframe_atg",
        }
        assert groups.map(lambda labels: labels == expected).all()
        assert grouped.size().eq(5).all()
        assert len(frame) == 10_000
    return len(frame)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--compare-reference", action="store_true")
    parser.add_argument("--skip-phylop", action="store_true")
    args = parser.parse_args()
    test_context_bounds()

    if args.compare_reference:
        compare_trees(
            args.root / "reference_data/regions",
            args.root / "Functional-Regions/regions",
        )
        for expected in sorted((args.root / "reference_data/atg").glob("chr*_atg_5way_labels.tsv")):
            actual = args.root / "ATG/source" / expected.name
            assert expected.read_bytes() == actual.read_bytes(), expected.name
        for name in ("all_chr_atg_5way.tsv", "sampled_examples_atg5.tsv"):
            assert (
                args.root / "reference_data/atg" / name
            ).read_bytes() == (
                args.root / "ATG/source" / name
            ).read_bytes(), name

    paths = [
        *[
            args.root / "Functional-Regions"
            / f"functional-{task}-gamba-{policy}.parquet"
            for task in ("upstream", "random", "random-noannot")
            for policy in ("causal", "bidi")
        ],
        *[
            args.root / "Functional-Regions"
            / f"functional-multiclass-gamba-{scope}-{policy}.parquet"
            for scope in ("full", "100bp")
            for policy in ("causal", "bidi")
        ],
        *[
            args.root / "ATG" / f"atg-gamba-{policy}.parquet"
            for policy in ("causal", "bidi")
        ],
    ]
    for path in paths:
        print(
            f"{path.relative_to(args.root)}: "
            f"{verify_parquet(path, not args.skip_phylop):,} rows"
        )
    print("all outputs verified")


def test_context_bounds() -> None:
    assert context_bounds(10_000, 5_000, 5_100, "+", 2_048, "causal") == (
        3_052,
        5_100,
    )
    assert context_bounds(10_000, 5_000, 5_100, "-", 2_048, "causal") == (
        5_000,
        7_048,
    )
    assert context_bounds(
        10_000, 5_000, 5_100, "+", 2_048, "symmetric"
    ) == (4_026, 6_074)


if __name__ == "__main__":
    main()
