#!/usr/bin/env python
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Train the XGBoost peak-filter used by ``isoquant_lib/terminal_counter.py``.

Inputs a per-peak feature CSV (the same format IsoQuant emits when training
data collection is enabled via the hidden ``--collect_polya_training`` /
``--collect_tss_training`` flags), trains an XGBoost classifier with the
hyperparameters used for the shipped model, and writes the JSON model file.

Typical use::

    # 1. Collect features from a real run (hidden dev-only flag):
    python isoquant.py --genedb ANNOTATION.gtf --bam READS.bam --data_type nanopore \\
        -o train_run --prefix run --collect_polya_training peaks.csv

    # 2. Train a fresh model:
    python misc/train_polya_tss_model.py \\
        --features peaks.csv \\
        --output isoquant_lib/data/model_polya.json

The features CSV must contain at least the columns referenced by
``FEATURE_COLUMNS`` plus a binary ``true_peak`` label column and a
``chromosome`` column for train/test splitting.
"""

import argparse
import logging
import random
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
from xgboost import XGBClassifier

# Importable so we always train on the same feature set the runtime expects.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from isoquant_lib.terminal_counter import FEATURE_COLUMNS  # noqa: E402

logger = logging.getLogger("train_polya_tss_model")

# Hyperparameters used for the shipped polyA / TSS models.
DEFAULT_PARAMS = {
    'booster': 'dart',
    'learning_rate': 0.25989024187450366,
    'gamma': 4.087502269121285,
    'max_depth': 6,
    'min_child_weight': 5.722951381412213,
    'reg_lambda': 4.84758559371017,
    'reg_alpha': 0.03795296107941127,
    'subsample': 0.8809954546952385,
    'colsample_bytree': 0.7483035972922923,
    'n_estimators': 111,
    'scale_pos_weight': 0.297061977156673,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True, type=Path,
                   help="CSV with per-peak features and 'true_peak' label.")
    p.add_argument("--output", required=True, type=Path,
                   help="Path to write the trained XGBoost JSON model.")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for the chromosome-level train/test split.")
    p.add_argument("--test-fraction", type=float, default=0.5,
                   help="Fraction of chromosomes to use for the test split.")
    return p.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()

    df = pd.read_csv(args.features)
    if 'true_peak' not in df.columns:
        raise SystemExit("features CSV is missing required 'true_peak' label column")
    missing = [c for c in FEATURE_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"features CSV is missing feature columns: {missing}")

    chrs = sorted(df['chromosome'].dropna().unique().tolist())
    if len(chrs) < 2:
        raise SystemExit("need at least two distinct chromosomes for a train/test split")
    rng = random.Random(args.seed)
    n_test = max(1, int(round(len(chrs) * args.test_fraction)))
    test_chrs = set(rng.sample(chrs, n_test))
    train_chrs = [c for c in chrs if c not in test_chrs]

    X = df[FEATURE_COLUMNS].astype(float)
    y = df['true_peak'].astype(int)
    train_mask = df['chromosome'].isin(train_chrs)
    test_mask = df['chromosome'].isin(test_chrs)

    model = XGBClassifier(**DEFAULT_PARAMS)
    model.fit(X[train_mask], y[train_mask])

    y_pred = model.predict(X[test_mask])
    y_score = model.predict_proba(X[test_mask])[:, 1]
    logger.info("test F1: %.4f", f1_score(y[test_mask], y_pred))
    try:
        logger.info("test ROC-AUC: %.4f", roc_auc_score(y[test_mask], y_score))
    except ValueError as e:
        logger.warning("ROC-AUC undefined: %s", e)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    model.save_model(str(args.output))
    logger.info("wrote model to %s", args.output)


if __name__ == "__main__":
    main()
