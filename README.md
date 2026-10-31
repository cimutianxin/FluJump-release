# FluJump-release

A minimal implementation of the influenza A hemagglutinin (HA) surveillance-scoring workflow described in the FluJump manuscript. It trains frozen ESM-2 linear readouts and converts their raw decision scores into training-set quantiles and alert tiers.

This repository contains source-data acquisition and the minimal modelling workflow, not the full research pipeline or experiment suite.

## Implemented workflow

- Raw, unaligned HA amino-acid sequences; standard 20-amino-acid alphabet and length 300–700 residues.
- Frozen `facebook/esm2_t30_150M_UR50D` representations, mean-pooled at layer 28 by default. Special tokens and padding are excluded.
- One ridge logistic-regression readout per binary label. `StandardScaler` is fitted on training data only.
- Five-fold stratified grid search over `C = {1e-3, 1e-2, 1e-1, 1, 10, 100}` with seed 42.
- Raw `decision_function` scores throughout; no sigmoid transformation or post-hoc sign flipping.
- Training quantiles and three tiers:
  - `high`: score ≥ median training-positive score.
  - `medium`: score ≥ 95th percentile of training-negative scores and below the high threshold.
  - `low`: all remaining scores.
- Scope flags: H1/H3/H5 are marked `validated`; H7/H10 are `outside_validated_scope`; other subtypes are `unvalidated`.

Scores and tiers are ranking aids, not transmission diagnoses or calibrated risk estimates.

## Scope

The release covers source-data acquisition, model fitting and scoring. It does **not** implement downstream data cleaning, CD-HIT clustering, cluster-label construction, target-side layer selection, subtype/year-forward evaluation, comparators, or site-attribution experiments.

Training labels must therefore be prepared beforehand. In the manuscript, `label_is_jump` denotes a cluster observed in at least two host or environmental surveillance categories, and `label_is_jump_human` denotes its subset with human plus at least one other non-human category. All isolates inherit their cluster label. Use only the intended H1+H3 training rows when reproducing the manuscript workflow; the script does not filter subtypes or splits automatically.

## Installation

Python 3.10 or newer is recommended.

```bash
git clone https://github.com/cimutianxin/FluJump-release.git
cd FluJump-release
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The first model run downloads the ESM-2 backbone from Hugging Face. If it is already cached, add `--offline` to either model command to disable network access.

## Source-data acquisition

Run all acquisition commands from the repository root. Downloaded files are written under `data/`, which is excluded from Git.

### NCBI Protein HA sequences and metadata

Set the email required by NCBI E-utilities, then run the original primary and H5 supplementary queries:

```bash
export NCBI_EMAIL="your-email@example.org"
python data_acquisition/download_flu.py
python data_acquisition/supplement_h5.py
```

`NCBI_API_KEY` can also be set when higher request limits are needed. The scripts write `data/raw/H1_raw.csv`, `H3_raw.csv`, `H5_raw.csv` and `H7_raw.csv`, including HA sequences, host annotations, collection dates, geography and other GenBank qualifiers. The H5 script adds bovine and human records and deduplicates by accession.

NCBI Protein is a live database, so a new download can differ from the manuscript collection. The fixed archived dataset, rather than a fresh query, is required for exact numerical reproduction.

### Borkenhagen receptor-binding data

Download Data S1 sequences and Data S2 labels from the open PMC archive of [Borkenhagen and Runstadler (2024)](https://doi.org/10.1111/irv.70044), then merge them using the identifier-matching procedure from the research repository:

```bash
python data_acquisition/download_borkenhagen.py
python data_acquisition/merge_borkenhagen.py
```

Raw files are saved in `data/dataset_borkenhagen_raw/`; the merged table is written to `data/dataset_borkenhagen/borkenhagen_merged.csv`. These commands acquire and merge the source data only; they do not recreate the manuscript's control filtering or stratified split.

## 1. Train

The training CSV needs one sequence column and at least one binary label column. The default column names are:

```text
sequence,label_is_jump,label_is_jump_human
```

Train both manuscript tasks and save one model bundle:

```bash
python flujump.py train \
  --train-csv training.csv \
  --model-out artifacts/flujump.joblib \
  --device auto
```

For differently named columns or a single task:

```bash
python flujump.py train \
  --train-csv training.csv \
  --model-out artifacts/cross_species.joblib \
  --sequence-column ha_sequence \
  --label-columns label_is_jump
```

Each class must contain at least five training rows. The selected regularization value and alert thresholds are printed during training. Model bundles can contain large scikit-learn objects; load only bundles from trusted sources.

## 2. Score a CSV

The input CSV must contain `sequence` and `subtype` columns by default. Subtypes may be written as `H5` or `H5N1`.

```bash
python flujump.py score \
  --model artifacts/flujump.joblib \
  --input-csv sequences.csv \
  --output-csv scores.csv \
  --device auto
```

The output preserves the input columns and adds, for every trained label:

- `raw_score__<label>`: untransformed linear-readout score.
- `train_quantile__<label>`: quantile in the training-score distribution, from 0 to 1.
- `tier__<label>`: `high`, `medium`, or `low`.

It also adds `tiering_scope` and `scope_note`. H7/H10 scores are still returned, but their tiers must not be interpreted as validated alerts.

## 3. Score one sequence

```bash
python flujump.py score \
  --model artifacts/flujump.joblib \
  --sequence "MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI" \
  --subtype H5N1
```

Add `--output-json score.json` to save the result instead of printing it.

## Reference environment

The manuscript analyses used Python 3.10, PyTorch 2.11, Hugging Face Transformers 5.13, and scikit-learn 1.7. This minimal repository exposes the core method rather than reproducing the manuscript's numerical results.
