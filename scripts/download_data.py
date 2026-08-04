#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import shutil
import urllib.request
from pathlib import Path


GAMBA_COMMIT = "e83984ea20bb4ee017993144fbe17e7bae3cdddc"
GAMBA_RAW = f"https://raw.githubusercontent.com/microsoft/gamba/{GAMBA_COMMIT}/data_processing/region_info"
GENCODE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
FASTA = "https://storage.googleapis.com/basenji_barnyard2/hg38.ml.fa.gz"
RMSK = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
REFGENE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
PHYLOP = "https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig"
PHYLOP_NAME = "241-mammalian-2020v2.bigWig"
PHYLOP_SIZE = 21_888_290_307
SMALL_FILES = [
    "experiments.tsv", "promoters.bed", "hg38_UCNE_coordinates.bed",
    "ucne_paralogues.txt",
]
EXPECTED_DERIVED = {
    "repeats_hg38.bed": "2e3f90cb003b45c81c25c5f29a94cc8b70908be0554e67c208ac8548996ca0c9",
    "UCSC_5UTR_exons.bed": "a6aeb685398856f186f83a4be4aa2a9745e40413079b87f6bd056ce7a9b097a2",
    "UCSC_3UTR_exons.bed": "31b56efb0e8e726318faeee3bb056578de228949131ab3b4d1d422d410ba1017",
}
EXPECTED_FASTA = "a3145abd8cad4cfdedde7e68a5c6263d4dfbdf1bfe4008beac9f002ddc6c8149"


def download(url: str, path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with urllib.request.urlopen(url) as source, temporary.open("wb") as target:
        shutil.copyfileobj(source, target)
    temporary.replace(path)


def verify(path: Path, expected: str) -> None:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected:
        raise RuntimeError(f"{path}: expected sha256 {expected}, got {actual}")


def verify_size(path: Path, expected: int) -> None:
    actual = path.stat().st_size
    if actual != expected:
        raise RuntimeError(f"{path}: expected {expected} bytes, got {actual}")


def gzip_lines(url: str):
    response = urllib.request.urlopen(url)
    return io.TextIOWrapper(gzip.GzipFile(fileobj=response))


def split_gtf(output: Path) -> None:
    if (output / "chr1.gtf").exists():
        return
    temporary = output.with_name(output.name + ".tmp")
    shutil.rmtree(temporary, ignore_errors=True)
    temporary.mkdir(parents=True)
    handles = {}
    try:
        with gzip_lines(GENCODE) as source:
            for line in source:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t", 1)[0]
                if not chrom.startswith("chr"):
                    continue
                if chrom not in handles:
                    handles[chrom] = (temporary / f"{chrom}.gtf").open("w")
                handles[chrom].write(line)
    finally:
        for handle in handles.values():
            handle.close()
    shutil.rmtree(output, ignore_errors=True)
    temporary.replace(output)


def build_repeats(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip_lines(RMSK) as source, temporary.open("w") as output:
        for fields in csv.reader(source, delimiter="\t"):
            output.write("\t".join([
                fields[5], fields[6], fields[7], fields[10], fields[1], fields[11],
            ]) + "\n")
    temporary.replace(path)


def build_utrs(output: Path) -> None:
    utr5, utr3 = output / "UCSC_5UTR_exons.bed", output / "UCSC_3UTR_exons.bed"
    if utr5.exists() and utr3.exists():
        return
    five_tmp, three_tmp = utr5.with_suffix(".bed.tmp"), utr3.with_suffix(".bed.tmp")
    with gzip_lines(REFGENE) as source, five_tmp.open("w") as five, three_tmp.open("w") as three:
        for fields in csv.reader(source, delimiter="\t"):
            name, chrom, strand = fields[1], fields[2], fields[3]
            tx_start, tx_end = int(fields[4]), int(fields[5])
            cds_start, cds_end = int(fields[6]), int(fields[7])
            exon_starts = list(map(int, fields[9].rstrip(",").split(",")))
            exon_ends = list(map(int, fields[10].rstrip(",").split(",")))
            ranges = (
                ((tx_start, cds_start), (cds_end, tx_end))
                if strand == "+" else ((cds_end, tx_end), (tx_start, cds_start))
            )
            for label, bounds, target in (
                ("utr5", ranges[0], five), ("utr3", ranges[1], three)
            ):
                for index, (start, end) in enumerate(zip(exon_starts, exon_ends)):
                    left, right = max(start, bounds[0]), min(end, bounds[1])
                    if left < right:
                        target.write(
                            f"{chrom}\t{left}\t{right}\t{name}_{label}_{index}\t0\t{strand}\n"
                        )
    five_tmp.replace(utr5)
    three_tmp.replace(utr3)


def symlink(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    target.unlink(missing_ok=True)
    target.symlink_to(source)


def reuse_existing(old_root: Path, data_root: Path, gamba_root: Path,
                   include_phylop: bool) -> None:
    old_data = old_root / "data_processing/data"
    old_info = old_root / "data_processing/region_info"
    symlink(old_data / "240-mammalian/hg38.ml.fa", data_root / "hg38.ml.fa")
    symlink(old_data / "240-mammalian/hg38.ml.fa.fai", data_root / "hg38.ml.fa.fai")
    symlink(old_data / "gtfs", data_root / "gtfs")
    if include_phylop:
        source = old_data / "240-mammalian" / PHYLOP_NAME
        if not source.exists():
            raise FileNotFoundError(source)
        symlink(source, data_root / PHYLOP_NAME)
    info = data_root / "region_info"
    for name in ("repeats_hg38.bed", "UCSC_5UTR_exons.bed", "UCSC_3UTR_exons.bed"):
        symlink(old_info / name, info / name)
    for name in SMALL_FILES:
        target = info / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((gamba_root / "data_processing/region_info" / name).read_bytes())


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=root / "data")
    parser.add_argument("--reuse-from", type=Path)
    parser.add_argument("--gamba-root", type=Path, default=root.parent / "Gamba")
    parser.add_argument(
        "--skip-phylop",
        action="store_true",
        help="Skip the default 21.9 GB Zoonomia phyloP bigWig.",
    )
    args = parser.parse_args()
    include_phylop = not args.skip_phylop

    if args.reuse_from:
        reuse_existing(
            args.reuse_from, args.data_root, args.gamba_root, include_phylop
        )
        verify(args.data_root / "hg38.ml.fa", EXPECTED_FASTA)
        for name, digest in EXPECTED_DERIVED.items():
            verify(args.data_root / "region_info" / name, digest)
        if include_phylop:
            verify_size(args.data_root / PHYLOP_NAME, PHYLOP_SIZE)
        print(f"linked existing inputs into {args.data_root}")
        return

    info = args.data_root / "region_info"
    for name in SMALL_FILES:
        download(f"{GAMBA_RAW}/{name}", info / name)

    fasta = args.data_root / "hg38.ml.fa"
    if not fasta.exists():
        compressed = args.data_root / "hg38.ml.fa.gz"
        download(FASTA, compressed)
        temporary = fasta.with_suffix(".fa.tmp")
        with gzip.open(compressed, "rb") as source, temporary.open("wb") as target:
            shutil.copyfileobj(source, target)
        temporary.replace(fasta)
        compressed.unlink()
    split_gtf(args.data_root / "gtfs")
    build_repeats(info / "repeats_hg38.bed")
    build_utrs(info)
    verify(fasta, EXPECTED_FASTA)
    for name, digest in EXPECTED_DERIVED.items():
        verify(info / name, digest)
    if include_phylop:
        download(PHYLOP, args.data_root / PHYLOP_NAME)
        verify_size(args.data_root / PHYLOP_NAME, PHYLOP_SIZE)
    print(f"downloaded and prepared inputs under {args.data_root}")


if __name__ == "__main__":
    main()
