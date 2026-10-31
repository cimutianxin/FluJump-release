#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

__version__ = "0.1.0"

DEFAULT_MODEL = "facebook/esm2_t30_150M_UR50D"
DEFAULT_LAYER = 28
DEFAULT_LABELS = ["label_is_jump", "label_is_jump_human"]
RIDGE_C_VALUES = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0]
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
VALIDATED_SUBTYPES = {"H1", "H3", "H5"}
OUTSIDE_SCOPE_SUBTYPES = {"H7", "H10"}


def read_csv(path: str):
    import pandas as pd

    csv_path = Path(path)
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    frame = pd.read_csv(csv_path)
    if frame.empty:
        raise ValueError(f"CSV file contains no rows: {csv_path}")
    return frame


def normalize_sequence(value, context: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{context}: sequence must be a string")
    sequence = value.strip().upper()
    if not sequence:
        raise ValueError(f"{context}: sequence is empty")
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if invalid:
        raise ValueError(
            f"{context}: non-standard amino-acid characters: {', '.join(invalid)}"
        )
    if not 300 <= len(sequence) <= 700:
        raise ValueError(
            f"{context}: expected an HA sequence of 300–700 residues, got {len(sequence)}"
        )
    return sequence


def get_sequences(frame, column: str) -> list[str]:
    if column not in frame.columns:
        raise ValueError(f"Missing sequence column: {column}")
    return [
        normalize_sequence(value, f"row {index + 2}")
        for index, value in enumerate(frame[column])
    ]


def get_labels(frame, column: str):
    import numpy as np
    import pandas as pd

    if column not in frame.columns:
        raise ValueError(f"Missing label column: {column}")
    if frame[column].isna().any():
        raise ValueError(f"Label column contains missing values: {column}")
    values = pd.to_numeric(frame[column], errors="raise").to_numpy(dtype=float)
    if not np.isin(values, [0.0, 1.0]).all():
        raise ValueError(f"Label column must contain only 0 and 1: {column}")
    labels = values.astype(int)
    counts = np.bincount(labels, minlength=2)
    if counts.min() < 5:
        raise ValueError(
            f"Label {column} needs at least five rows per class for five-fold CV; "
            f"got negative={counts[0]}, positive={counts[1]}"
        )
    return labels


def resolve_device(requested: str):
    import torch

    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


class ESM2Embedder:
    def __init__(self, model_name: str, layer: int, device: str):
        import torch
        from transformers import AutoTokenizer, EsmModel

        self.torch = torch
        self.device = resolve_device(device)
        self.layer = layer
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = EsmModel.from_pretrained(
            model_name, output_hidden_states=True
        ).to(self.device).eval()
        layer_count = int(self.model.config.num_hidden_layers)
        if not 0 <= layer <= layer_count:
            raise ValueError(
                f"Layer must be between 0 and {layer_count} for {model_name}; got {layer}"
            )
        print(f"Loaded {model_name} on {self.device}; layer={layer}", file=sys.stderr)

    def embed(self, sequences: list[str], batch_size: int):
        import numpy as np

        if batch_size < 1:
            raise ValueError("Batch size must be at least 1")
        batches = []
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]
            inputs = self.tokenizer(
                batch,
                padding=True,
                return_tensors="pt",
                return_special_tokens_mask=True,
            )
            model_inputs = {
                "input_ids": inputs["input_ids"].to(self.device),
                "attention_mask": inputs["attention_mask"].to(self.device),
            }
            residue_mask = (
                inputs["attention_mask"].bool()
                & ~inputs["special_tokens_mask"].bool()
            ).to(self.device)
            with self.torch.inference_mode():
                outputs = self.model(**model_inputs)
                hidden = outputs.hidden_states[self.layer]
                lengths = residue_mask.sum(dim=1)
                pooled = (hidden * residue_mask.unsqueeze(-1)).sum(dim=1)
                pooled = pooled / lengths.unsqueeze(-1)
            batches.append(pooled.cpu().float().numpy())
            print(
                f"Embedded {min(start + batch_size, len(sequences))}/{len(sequences)}",
                file=sys.stderr,
                flush=True,
            )
        return np.concatenate(batches, axis=0).astype(np.float32, copy=False)


def fit_readout(features, labels):
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GridSearchCV, StratifiedKFold
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    scaled = scaler.fit_transform(features)
    search = GridSearchCV(
        LogisticRegression(
            penalty="l2", solver="lbfgs", max_iter=5000, random_state=42
        ),
        {"C": RIDGE_C_VALUES},
        cv=StratifiedKFold(5, shuffle=True, random_state=42),
        scoring="roc_auc",
        n_jobs=-1,
    )
    search.fit(scaled, labels)
    classifier = search.best_estimator_
    training_scores = classifier.decision_function(scaled)
    high = float(np.median(training_scores[labels == 1]))
    medium = float(np.quantile(training_scores[labels == 0], 0.95))
    if high <= medium:
        raise ValueError(
            f"Degenerate alert thresholds: high={high:.6g}, medium={medium:.6g}"
        )
    return {
        "scaler": scaler,
        "classifier": classifier,
        "best_c": float(search.best_params_["C"]),
        "training_scores": np.sort(training_scores),
        "high_threshold": high,
        "medium_threshold": medium,
    }


def train_command(args) -> None:
    import joblib

    frame = read_csv(args.train_csv)
    sequences = get_sequences(frame, args.sequence_column)
    labels = {name: get_labels(frame, name) for name in args.label_columns}
    embedder = ESM2Embedder(args.model_name, args.layer, args.device)
    features = embedder.embed(sequences, args.batch_size)
    tasks = {}
    for name, task_labels in labels.items():
        task = fit_readout(features, task_labels)
        tasks[name] = task
        print(
            f"{name}: n={len(task_labels)}, positives={int(task_labels.sum())}, "
            f"best_C={task['best_c']:g}, high={task['high_threshold']:.6g}, "
            f"medium={task['medium_threshold']:.6g}",
            file=sys.stderr,
        )
    bundle = {
        "format_version": 1,
        "model_name": args.model_name,
        "layer": args.layer,
        "tasks": tasks,
    }
    output = Path(args.model_out)
    output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, output)
    print(f"Saved model bundle: {output}", file=sys.stderr)


def subtype_scope(subtype: str) -> dict[str, str]:
    normalized = str(subtype).strip().upper()
    match = re.match(r"^(H\d+)", normalized)
    ha_subtype = match.group(1) if match else normalized
    if ha_subtype in OUTSIDE_SCOPE_SUBTYPES:
        return {
            "ha_subtype": ha_subtype,
            "status": "outside_validated_scope",
            "note": "H7/H10 tiers are outside the validated application scope.",
        }
    if ha_subtype in VALIDATED_SUBTYPES:
        return {
            "ha_subtype": ha_subtype,
            "status": "validated",
            "note": "",
        }
    return {
        "ha_subtype": ha_subtype,
        "status": "unvalidated",
        "note": "This subtype has not been validated for alert tiering.",
    }


def score_embedding(task: dict, embedding) -> dict[str, float | str]:
    import numpy as np

    score = float(
        task["classifier"].decision_function(
            task["scaler"].transform(embedding.reshape(1, -1))
        )[0]
    )
    training_scores = task["training_scores"]
    quantile = float(
        np.searchsorted(training_scores, score, side="left") / len(training_scores)
    )
    if score >= task["high_threshold"]:
        tier = "high"
    elif score >= task["medium_threshold"]:
        tier = "medium"
    else:
        tier = "low"
    return {"raw_score": score, "train_quantile": quantile, "tier": tier}


def load_bundle(path: str) -> dict:
    import joblib

    model_path = Path(path)
    if not model_path.is_file():
        raise FileNotFoundError(f"Model bundle not found: {model_path}")
    bundle = joblib.load(model_path)
    required = {"format_version", "model_name", "layer", "tasks"}
    if not isinstance(bundle, dict) or not required.issubset(bundle):
        raise ValueError(f"Invalid FluJump model bundle: {model_path}")
    if bundle["format_version"] != 1 or not bundle["tasks"]:
        raise ValueError(f"Unsupported or empty FluJump model bundle: {model_path}")
    return bundle


def score_command(args) -> None:
    bundle = load_bundle(args.model)
    if args.sequence is not None:
        sequence = normalize_sequence(args.sequence, "command-line sequence")
        embedder = ESM2Embedder(bundle["model_name"], bundle["layer"], args.device)
        embedding = embedder.embed([sequence], 1)[0]
        result = {
            "subtype": args.subtype,
            "scope": subtype_scope(args.subtype),
            "model_name": bundle["model_name"],
            "layer": bundle["layer"],
            "tasks": {
                name: score_embedding(task, embedding)
                for name, task in bundle["tasks"].items()
            },
        }
        rendered = json.dumps(result, indent=2, ensure_ascii=False)
        if args.output_json:
            output = Path(args.output_json)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(rendered + "\n", encoding="utf-8")
            print(f"Saved scores: {output}", file=sys.stderr)
        else:
            print(rendered)
        return

    frame = read_csv(args.input_csv)
    sequences = get_sequences(frame, args.sequence_column)
    if args.subtype_column not in frame.columns:
        raise ValueError(f"Missing subtype column: {args.subtype_column}")
    if frame[args.subtype_column].isna().any():
        raise ValueError(f"Subtype column contains missing values: {args.subtype_column}")
    subtypes = frame[args.subtype_column].astype(str).tolist()
    embedder = ESM2Embedder(bundle["model_name"], bundle["layer"], args.device)
    features = embedder.embed(sequences, args.batch_size)
    output_frame = frame.copy()
    scopes = [subtype_scope(subtype) for subtype in subtypes]
    output_frame["tiering_scope"] = [scope["status"] for scope in scopes]
    output_frame["scope_note"] = [scope["note"] for scope in scopes]
    for name, task in bundle["tasks"].items():
        scores = [score_embedding(task, embedding) for embedding in features]
        output_frame[f"raw_score__{name}"] = [score["raw_score"] for score in scores]
        output_frame[f"train_quantile__{name}"] = [
            score["train_quantile"] for score in scores
        ]
        output_frame[f"tier__{name}"] = [score["tier"] for score in scores]
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    output_frame.to_csv(output, index=False)
    print(f"Saved {len(output_frame)} scored rows: {output}", file=sys.stderr)


def add_runtime_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda, or mps")
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Load the ESM-2 backbone only from the local Hugging Face cache",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Minimal frozen ESM-2 linear-readout workflow for influenza HA"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    train = subparsers.add_parser("train", help="Train linear readouts from a CSV")
    train.add_argument("--train-csv", required=True)
    train.add_argument("--model-out", required=True, help="Output .joblib bundle")
    train.add_argument("--sequence-column", default="sequence")
    train.add_argument("--label-columns", nargs="+", default=DEFAULT_LABELS)
    train.add_argument("--model-name", default=DEFAULT_MODEL)
    train.add_argument("--layer", type=int, default=DEFAULT_LAYER)
    add_runtime_arguments(train)
    train.set_defaults(handler=train_command)

    score = subparsers.add_parser("score", help="Score one sequence or a CSV")
    score.add_argument("--model", required=True, help="Trained .joblib bundle")
    source = score.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-csv", help="CSV containing sequences and subtypes")
    source.add_argument("--sequence", help="One raw, unaligned HA amino-acid sequence")
    score.add_argument("--subtype", help="Subtype for --sequence, for example H5N1")
    score.add_argument("--sequence-column", default="sequence")
    score.add_argument("--subtype-column", default="subtype")
    score.add_argument("--output-csv", help="Required with --input-csv")
    score.add_argument("--output-json", help="Optional with --sequence")
    add_runtime_arguments(score)
    score.set_defaults(handler=score_command)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "score":
        if args.input_csv:
            if not args.output_csv:
                parser.error("--output-csv is required with --input-csv")
            if args.output_json:
                parser.error("--output-json is only used with --sequence")
        if args.sequence is not None:
            if not args.subtype:
                parser.error("--subtype is required with --sequence")
            if args.output_csv:
                parser.error("--output-csv is only used with --input-csv")
    if getattr(args, "offline", False):
        os.environ["HF_HUB_OFFLINE"] = "1"
        os.environ["TRANSFORMERS_OFFLINE"] = "1"
    try:
        args.handler(args)
    except (FileNotFoundError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")


if __name__ == "__main__":
    main()
