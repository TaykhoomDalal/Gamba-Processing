from __future__ import annotations

import hashlib
from collections.abc import Iterable
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SCHEMA = pa.schema([
    ("split", pa.string()),
    ("sequence", pa.string()),
    ("label", pa.string()),
    ("pair_id", pa.string()),
    ("category", pa.string()),
    ("scope", pa.string()),
    ("chrom", pa.string()),
    ("start", pa.int64()),
    ("end", pa.int64()),
    ("strand", pa.string()),
    ("context_start", pa.int64()),
    ("context_end", pa.int64()),
    ("roi_start", pa.int32()),
    ("roi_end", pa.int32()),
    ("pool_start", pa.int32()),
    ("pool_end", pa.int32()),
    ("name", pa.string()),
])
ATG_SCHEMA = pa.schema([
    *SCHEMA,
    ("transcript_id", pa.string()),
    ("gene_id", pa.string()),
    ("label_id", pa.int8()),
    ("delta_bp", pa.int64()),
])
PHYLOP_NAMES = (
    "mean", "std", "frac_pos", "frac_neg", "mean_pos", "mean_neg"
)
PHYLOP_FIELDS = tuple(
    pa.field(f"{prefix}_{name}", pa.float32())
    for prefix in ("phylop", "phylop_context")
    for name in PHYLOP_NAMES
)


def with_phylop(schema: pa.Schema) -> pa.Schema:
    return pa.schema([*schema, *PHYLOP_FIELDS])


def context_for_region(genome, chrom: str, start: int, end: int, strand: str,
                       length: int = 2048) -> dict | None:
    start, end = sorted((int(start), int(end)))
    feature_length = end - start
    if feature_length <= 0 or chrom not in genome:
        return None

    chrom_length = len(genome[chrom])
    if feature_length > length:
        kept = min(1000, length)
        if strand == "+":
            context_end = end
            context_start = max(0, context_end - kept)
        else:
            context_start, context_end = start, min(chrom_length, start + kept)
    elif strand != "+":
        context_start = start
        context_end = min(chrom_length, start + length)
    else:
        context_end = end
        context_start = max(0, end - length)

    sequence = genome[chrom][context_start:context_end].seq.upper()
    roi_start = max(0, start - context_start)
    roi_end = min(len(sequence), end - context_start)
    if roi_end <= roi_start:
        return None
    if strand == "-":
        sequence = sequence.translate(COMPLEMENT)[::-1]
        roi_start, roi_end = len(sequence) - roi_end, len(sequence) - roi_start

    return {
        "sequence": sequence,
        "context_start": context_start,
        "context_end": context_end,
        "roi_start": roi_start,
        "roi_end": roi_end,
    }


def stable_100bp_span(category: str, pair_id: str, start: int, end: int,
                      seed: int) -> tuple[int, int] | None:
    if end - start < 100:
        return None
    key = f"{seed}|{category}|{pair_id}".encode()
    offset = int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "little")
    span_start = start + offset % (end - start - 99)
    return span_start, span_start + 100


def write_parquet(path: Path, rows: Iterable[dict], schema: pa.Schema = SCHEMA,
                  batch_size: int = 2048) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    writer = pq.ParquetWriter(temporary, schema, compression="zstd")
    batch, count = [], 0
    try:
        for row in rows:
            batch.append({field.name: row.get(field.name) for field in schema})
            if len(batch) == batch_size:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                count += len(batch)
                batch.clear()
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
            count += len(batch)
    finally:
        writer.close()
    temporary.replace(path)
    return count
