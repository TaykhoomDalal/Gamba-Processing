# GAMBA evaluation parquet processing

This repository is the small, data-only part of
[Microsoft GAMBA](https://github.com/microsoft/gamba): it recreates the
functional-region and ATG evaluation inputs and writes fixed, model-independent
parquets. It contains no model, embedding, plotting, or metric code.

Each logical dataset is emitted twice:

- `*-causal.parquet`: ROI end-anchored after strand orientation for
  autoregressive models such as Evo2;
- `*-bidi.parquet`: ROI centered for bidirectional models such as the
  distilled student, GPN-Star, and PhyloGPN.

The implementation is derived from GAMBA commit
[`e83984e`](https://github.com/microsoft/gamba/commit/e83984ea20bb4ee017993144fbe17e7bae3cdddc),
principally:

- `data_processing/create_eval_data.py`
- `data_processing/make_ATG_data.py`
- `src/evaluation/run_eval.py`
- `src/evaluation/ATG_reps.py`
- `src/evaluation/utils/helpers.py`

GAMBA is MIT licensed; the license and attribution are retained in
[LICENSE](LICENSE) and [NOTICE.md](NOTICE.md).

## Outputs

```text
Functional-Regions/
├── functional-upstream-gamba-{causal,bidi}.parquet
├── functional-random-gamba-{causal,bidi}.parquet
├── functional-random-noannot-gamba-{causal,bidi}.parquet
├── functional-multiclass-gamba-full-{causal,bidi}.parquet
├── functional-multiclass-gamba-100bp-{causal,bidi}.parquet
├── regions/
└── regions-noncoding-added/

ATG/
├── atg-gamba-{causal,bidi}.parquet
└── source/
    ├── chr1_atg_5way_labels.tsv
    ├── ...
    ├── chr22_atg_5way_labels.tsv
    ├── all_chr_atg_5way.tsv
    └── sampled_examples_atg5.tsv
```

The functional parquets include `chr1`–`chr22` and `chrX`. Rows on `chr2`,
`chr3`, `chr16`, and `chr22` are marked `test` because those chromosomes were
held out from GAMBA pretraining. The remaining chromosomes are marked `train`.
Filtering to `split == "test"` reproduces the chromosome subset used by
GAMBA's representation evaluator.

The Hugging Face configs expose each physical parquet as the Hub split `all`.
That is separate from the parquet's `split` column: `all` means "load this
whole file," while the column contains the chromosome-level `train`/`test`
assignment. This avoids storing duplicate physical train and test files.

| Parquet | Labels |
|---|---|
| `functional-upstream-gamba*.parquet` | `feature`, `upstream` |
| `functional-random-gamba*.parquet` | `feature`, `random` |
| `functional-random-noannot-gamba*.parquet` | `feature`, `random-noannot` |
| `functional-multiclass-gamba-{full,100bp}*.parquet` | one label per functional category |
| `ATG/atg-gamba-{causal,bidi}.parquet` | `start`, `noncoding_near`, `noncoding_far`, `inframe_methionine`, `outframe_atg` |

The multiclass full-region and 100 bp variants are separate parquet files and
separate Hugging Face configs. Every functional parquet is a strict superset:
it includes the explicit `noncoding_regions` extension, while
`category != "noncoding_regions"` is exactly the canonical ten-category GAMBA
paper dataset. This avoids storing the canonical rows twice.

## Parquet schemas

All functional parquets use exactly the same columns and Arrow types:

| Column | Type | Meaning |
|---|---|---|
| `split`, `sequence`, `label` | string | probe/fine-tuning `train` or held-out `test`; model input and target |
| `pair_id` | string | matched functional pair or stable feature identifier |
| `category`, `scope` | string | functional class and pooling scope |
| `context_policy` | string | `causal` or `symmetric` window placement |
| `chrom`, `start`, `end`, `strand` | string/int64 | 0-based, half-open feature coordinates |
| `context_start`, `context_end` | int64 | forward-genome coordinates of `sequence` |
| `roi_start`, `roi_end` | int32 | feature offsets in the strand-oriented sequence |
| `pool_start`, `pool_end` | int32 | offsets that the evaluator should pool |
| `name` | string | source annotation identifier |
| `phylop_{mean,std,frac_pos,frac_neg,mean_pos,mean_neg}` | float32 | GAMBA phyloP baseline over the task pooling span |
| `phylop_context_{mean,std,frac_pos,frac_neg,mean_pos,mean_neg}` | float32 | GAMBA phyloP baseline over its symmetric 2,048 bp context |

The ATG parquet uses the same core columns and adds `transcript_id`, `gene_id`,
`label_id`, and `delta_bp`. Functional parquets do not contain empty ATG-only
columns.

ATG stores five source contexts per transcript. For the paper leaderboard,
merge `noncoding_near` and `noncoding_far` into one `noncoding` class and use
cosine leave-one-out 1-nearest-neighbor balanced accuracy. The author
extraction script's five-way Euclidean output is a separate diagnostic.

The GAMBA paper uses model-appropriate context geometry. `-causal` sequences
place the ROI at the end of the strand-oriented 2,048 bp window so an
autoregressive model can condition on the preceding bases; use these for Evo2
and other left-to-right models. `-bidi` sequences center the ROI so both
flanks are visible; use these for GAMBA encoders, the distilled student,
GPN-Star, PhyloGPN, and other masked or bidirectional models. The
`context_policy` column is `causal` or `symmetric`, respectively.

Only `-` records are reverse complemented. Causal features longer than 2,048
bp retain GAMBA's 1,000 bp truncation; bidi files retain the centered 2,048 bp
span.

## Environment

```bash
bash scripts/create_environment.sh
mamba activate gamba-processing
```

The environment pins Python 3.12.13, NumPy 2.5.1, pandas 3.0.5, PyArrow
25.0.0, pyfaidx 0.9.0.4, and pyBigWig 0.3.25.

## Inputs

For a fresh standalone build:

```bash
python scripts/download_data.py
```

This downloads or derives:

| Input | Source |
|---|---|
| bundled VISTA, promoter, UCNE, and paralogue files | `https://raw.githubusercontent.com/microsoft/gamba/e83984ea20bb4ee017993144fbe17e7bae3cdddc/data_processing/region_info/<file>` |
| `hg38.ml.fa` | `https://storage.googleapis.com/basenji_barnyard2/hg38.ml.fa.gz` |
| per-chromosome GTFs | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz` |
| RepeatMasker BED | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz` |
| 5′/3′ UTR BEDs | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz` |
| 241-mammalian phyloP bigWig | `https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig` |

The RepeatMasker and UTR conversions reproduce the input layout used for the
local GAMBA run. Their expected SHA-256 hashes are checked after download so a
mutable UCSC table cannot silently change the benchmark. Downloaded and
generated data are gitignored.

To reuse the already downloaded files from the former evaluation workspace:

```bash
python scripts/download_data.py \
  --reuse-from /home/t-tdalal/glm/model/Gamba-evals-old \
  --gamba-root /home/t-tdalal/glm/model/Gamba
```

This only creates local symlinks for the large files.

## GAMBA phyloP baseline features

GAMBA's phyloP baseline converts per-base Zoonomia scores into six float32
features. All generated parquets include these features by default.

| Column suffix | Calculation |
|---|---|
| `mean` | mean phyloP score |
| `std` | population standard deviation |
| `frac_pos` | fraction of bases with score greater than zero |
| `frac_neg` | fraction of bases with score less than zero |
| `mean_pos` | mean over positive scores, or zero when absent |
| `mean_neg` | mean over negative scores, or zero when absent |

The implementation matches GAMBA's behavior: uncovered bigWig positions become
zero, scores are rounded to two decimal places before summarization, and the
statistics are stored as float32. Minus-strand score arrays are reversed before
reduction, matching GAMBA's operation order exactly.

The default data setup downloads the 21.9 GB bigWig and verifies its expected
21,888,290,307-byte size:

```bash
python scripts/download_data.py
```

The reuse workflow links the existing local copy by default:

```bash
python scripts/download_data.py \
  --reuse-from /home/t-tdalal/glm/model/Gamba-evals-old \
  --gamba-root /home/t-tdalal/glm/model/Gamba
```

The normal functional and ATG processors append both feature sets directly to
their standard parquet filenames. To create sequence-only parquets instead:

```bash
bash scripts/run_all.sh --skip-phylop
```

The default parquets contain:

- `phylop_{mean,std,frac_pos,frac_neg,mean_pos,mean_neg}` for GAMBA's baseline
  ROI. This is the annotated ROI for ordinary features, the centered symmetric
  2,048 bp slice for features at least 2,048 bp long, the fixed pooling span for
  100 bp configs, and the three-base ATG.
- `phylop_context_{mean,std,frac_pos,frac_neg,mean_pos,mean_neg}` for GAMBA's
  symmetric 2,048 bp phyloP baseline context.

`scripts/add_phylop.py` is also available for existing sequence-only parquets.
Use `--in-place` to replace them, or omit it to write
`phylop/<input-name>-phylop.parquet`.

## Build

Run both processors:

```bash
bash scripts/run_all.sh
```

Or independently:

```bash
python Functional-Regions/process.py \
  --include-noncoding \
  --regions-dir Functional-Regions/regions-noncoding-added
python ATG/process.py
```

Each processor includes phyloP by default. Pass `--skip-phylop` to either
processor, or to `scripts/run_all.sh`, for sequence-only outputs.

Functional regions are capped at 10,000 retained anchors per category with
seed 42. ATG generation follows GAMBA's MANE Select CDS logic, then samples
2,000 complete examples approximately evenly across `chr1`–`chr22`, also with
seed 42. Omit `--include-noncoding` only when regenerating the ten-category
reference BED set; published parquets use the noncoding-inclusive superset.

### Functional chromosome split

GAMBA's original pretraining split was:

- training: `chr1`, `chr4`–`chr15`, `chr17`–`chr21`, and `chrX`;
- validation: `chr3`, `chr16`;
- test: `chr2`, `chr22`.

`src/evaluation/run_eval.py` combines the validation and test chromosomes into
one four-chromosome evaluation group. The parquets therefore label all four as
`test` and label the remaining chromosomes as `train`.

For zero-shot evaluation, use `test` to reproduce the GAMBA paper or evaluate
all rows for broader chromosome coverage. For probing or supervised
fine-tuning, `train` can fit the task-specific head and `test` can evaluate it
on held-out chromosomes. The split also records GAMBA's pretraining exposure:
`train` chromosomes were seen during GAMBA pretraining, while `test`
chromosomes were not. In the paper's zero-shot evaluation, the GAMBA models
were frozen and evaluated on this held-out `test` subset without fitting a
task-specific model on these benchmark rows. The ATG benchmark follows GAMBA's
separate all-autosome sampling and remains a single `test` split.

## Validation

```bash
python scripts/verify_outputs.py
```

The completed reference comparison found:

- all 920 ten-category BED files byte-identical to unmodified GAMBA output;
- all 22 ATG chromosome TSVs, the combined TSV, and the sampled 2,000 examples
  byte-identical to GAMBA output;
- 169,674 paper rows in each binary parquet after excluding
  `noncoding_regions`;
- 188,562 total rows in each binary superset (`150,728` train and `37,834`
  test);
- 84,837 paper full-region multiclass rows and 58,643 paper 100 bp rows after
  the same filter;
- 94,281 total full-region rows and 67,624 total 100 bp rows;
- 10,000 rows in each ATG context-policy parquet.

## Intentional, visible differences from GAMBA

1. GAMBA's checked-in `experiments.tsv` uses `Element Coordinates`, while
   `create_eval_data.py` only accepts `coordinate_hg38` or `coord`. This parser
   accepts all three names; coordinates are transformed identically.
2. `create_eval_data.py` generates ten functional categories, but
   `run_eval.py` still lists `noncoding_regions`. Passing
   `--include-noncoding` adds a deterministic, half-open noncoding complement;
   without that optional flag the output is GAMBA's original ten categories.
3. GAMBA selects `roi100bp` with Python's process-randomized `hash()`. Here the
   same seed/category/pair ID is passed through BLAKE2 so the saved eval set is
   stable across machines and Python processes.
4. GAMBA stores a minus-strand ATG position as `[position, position + 3)` and
   then reverse complements that interval. This pipeline intentionally keeps
   that convention even when the resulting displayed triplet is surprising.

These are processing fixes, not model changes. The ten-category reference BED
output is used for direct equivalence checks against the unmodified GAMBA
generator; `scripts/run_all.sh` appends the noncoding category without
changing any canonical row.
