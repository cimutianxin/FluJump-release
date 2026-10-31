# FluJump-release

A minimal implementation of the influenza A hemagglutinin (HA) surveillance-scoring workflow described in the FluJump manuscript. It trains frozen ESM-2 linear readouts and converts their raw decision scores into training-set quantiles and alert tiers.

This repository intentionally contains no research-data pipeline or experiment suite. For the full analyses, see [cimutianxin/FluJump](https://github.com/cimutianxin/FluJump).

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

The minimal implementation covers model fitting and scoring only. It does **not** implement NCBI downloads, metadata standardization, CD-HIT clustering, cluster-label construction, target-side layer selection, subtype/year-forward evaluation, comparators, or site-attribution experiments.

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

The first run downloads the ESM-2 backbone from Hugging Face. If it is already cached, add `--offline` to either command to disable network access.

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
