# GAMBA evaluation parquet processing

This repository is the small, data-only part of
[Microsoft GAMBA](https://github.com/microsoft/gamba): it recreates the
functional-region and ATG evaluation inputs and writes fixed, model-independent
parquets. It contains no model, embedding, plotting, or metric code.

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
├── functional-upstream-gamba.parquet
├── functional-upstream-gamba-noncoding-added.parquet
├── functional-random-gamba.parquet
├── functional-random-gamba-noncoding-added.parquet
├── functional-random-noannot-gamba.parquet
├── functional-random-noannot-gamba-noncoding-added.parquet
├── functional-multiclass-gamba-full.parquet
├── functional-multiclass-gamba-100bp.parquet
├── functional-multiclass-gamba-full-noncoding-added.parquet
├── functional-multiclass-gamba-100bp-noncoding-added.parquet
├── regions/
└── regions-noncoding-added/

ATG/
├── atg-gamba.parquet
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
| `ATG/atg-gamba.parquet` | `start`, `noncoding_near`, `noncoding_far`, `inframe_methionine`, `outframe_atg` |

The multiclass full-region and 100 bp variants are separate parquet files and
separate Hugging Face configs. Files ending in `-noncoding-added.parquet`
include the explicit `noncoding_regions` extension; unsuffixed files are the
canonical ten-category GAMBA data.

## Parquet schemas

All ten functional parquets use exactly the same columns and Arrow types:

| Column | Type | Meaning |
|---|---|---|
| `split`, `sequence`, `label` | string | probe/fine-tuning `train` or held-out `test`; model input and target |
| `pair_id` | string | matched functional pair or stable feature identifier |
| `category`, `scope` | string | functional class and pooling scope |
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

Sequences use GAMBA's asymmetric 2,048 bp window: the ROI is at the right edge
on `+` records and at the left edge on every other strand value. Only `-`
records are reverse complemented. Features longer than 2,048 bp use GAMBA's
1,000 bp truncation.

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
python Functional-Regions/process.py
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
seed 42. Omit `--include-noncoding` to build the exact ten-category GAMBA
region set instead of the final extended parquet set.

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
- 169,674 rows in each canonical binary parquet
  (`135,692` train, `33,982` test);
- 188,560 rows in each binary parquet with noncoding added
  (`150,726` train, `37,834` test);
- 84,837 canonical full-region multiclass rows and 58,643 canonical 100 bp rows;
- 94,280 full-region and 67,623 100 bp rows with noncoding added;
- 10,000 rows in the ATG parquet.

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

These are processing fixes, not model changes. The ten-category
default BED output is used for direct equivalence checks against the
unmodified GAMBA generator; `scripts/run_all.sh` builds both canonical and
`noncoding-added` parquet variants.
