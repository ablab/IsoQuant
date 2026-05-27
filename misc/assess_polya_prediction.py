#!/usr/bin/env python3
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Evaluate IsoQuant's polyA-site predictions against a ground-truth GTF.

Designed to be called by ``isoquant_tests/github/run_pipeline.py`` under the
``polya_prediction`` run type. Emits a TSV with one ``metric<TAB>value`` row
per metric so ``check_value()`` in the launcher can compare against the
baselines block in the YAML config.

Ground truth comes from a GTF (typically a synthetic GTF with alternative
transcript ends used to simulate the reads). The GTF is loaded via
``gffutils`` and a ``.db`` cache is built next to it on first use so the
~156 K-record file is only parsed once per CI host. Prediction file is the
six-column ungrouped TSV emitted by ``isoquant_lib/terminal_counter.py``.

Annotated polyA position formula matches
``isoquant_lib.terminal_counter.PolyACounter._annotated_position``:
``last_exon_end + 1`` on ``+`` strand, ``first_exon_start - 1`` on ``-``.
"""

import argparse
import logging
import os
import sys
from collections import defaultdict
from traceback import print_exc
from typing import Optional

import gffutils

logger = logging.getLogger("assess_polya_prediction")

# IsoQuant's prediction TSV column layout (see EMPTY_COLUMNS in
# isoquant_lib/terminal_counter.py). The grouped variant adds two more
# columns; this script ignores them.
TSV_COLUMNS = ["chromosome", "transcript_id", "gene_id",
               "prediction", "counts", "flag"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prediction", required=True,
                   help="IsoQuant polyA_prediction.tsv (or .gz)")
    p.add_argument("--reference_gtf", required=True,
                   help="Ground-truth GTF (synthetic) or .db file")
    p.add_argument("--output", "-o", required=True,
                   help="Output TSV with one metric<TAB>value per row")
    p.add_argument("--transcript_counts",
                   help="IsoQuant transcript_counts.tsv — when provided, "
                        "only transcripts with count > 0 are included in "
                        "the recall denominator (filters out unexpressed).")
    p.add_argument("--tolerance", type=int, default=10,
                   help="bp tolerance for matching a prediction to a truth "
                        "end coordinate. Default 10 matches "
                        "ANNOTATION_TOLERANCE in terminal_counter.py.")
    return p.parse_args()


def get_or_build_db(gtf_path: str) -> gffutils.FeatureDB:
    """Return a FeatureDB for the supplied GTF, building a sibling .db cache
    on first call. Mirrors the flags used by ``isoquant_lib/gtf2db.py`` so
    the cache stays compatible with any other consumer in the project."""
    if gtf_path.endswith(".db"):
        return gffutils.FeatureDB(gtf_path, keep_order=True)
    db_path = gtf_path + ".db"
    if not os.path.exists(db_path):
        logger.info("Building gffutils DB at %s (one-time cache)", db_path)
        gffutils.create_db(
            gtf_path, db_path,
            force=False, keep_order=True, merge_strategy="error",
            sort_attribute_values=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
        )
    return gffutils.FeatureDB(db_path, keep_order=True)


def annotated_polya(strand: str, exon_blocks: list) -> Optional[int]:
    """Match PolyACounter._annotated_position. exon_blocks is a list of
    ``(start, end)`` tuples sorted by genomic coordinate."""
    if not exon_blocks:
        return None
    if strand == "+":
        return exon_blocks[-1][1] + 1
    if strand == "-":
        return exon_blocks[0][0] - 1
    return None


def _is_novel_id(tid: str) -> bool:
    """True when *tid* looks like ``ENSMUST00000000033.12_142204503``
    (base GENCODE ID + ``_`` + novel polyA coordinate)."""
    parts = tid.rsplit("_", 1)
    return len(parts) == 2 and parts[1].isdigit()


def _base_id(tid: str) -> str:
    return tid.rsplit("_", 1)[0]


def load_truth(db: gffutils.FeatureDB) -> dict:
    """Return ``{transcript_id: [(chr, polya_pos, is_novel), ...]}`` keyed
    by base GENCODE ID. Novel transcript entries (IDs with a ``_<digits>``
    suffix) are merged into the base transcript's list so that a prediction
    for the base ID can match either the known or the novel site."""
    truth: dict = defaultdict(list)
    n_known = n_novel = 0
    for t in db.features_of_type("transcript"):
        tid = t.id
        if not tid:
            continue
        exon_blocks = sorted(
            ((e.start, e.end) for e in db.children(t, featuretype="exon")),
            key=lambda b: b[0])
        polya = annotated_polya(t.strand, exon_blocks)
        if polya is None:
            continue
        novel = _is_novel_id(tid)
        key = _base_id(tid) if novel else tid
        truth[key].append((t.seqid, polya, novel))
        if novel:
            n_novel += 1
        else:
            n_known += 1
    logger.info("Loaded %d known + %d novel truth entries (%d base IDs)",
                n_known, n_novel, len(truth))
    return truth


def load_expressed_transcripts(counts_path: str) -> set:
    """Return the set of *base* transcript IDs with count > 0. Novel IDs
    (``base_<digits>``) are mapped to their base GENCODE ID."""
    expressed: set = set()
    with open(counts_path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and float(parts[1]) > 0:
                tid = parts[0]
                expressed.add(_base_id(tid) if _is_novel_id(tid) else tid)
    logger.info("Loaded %d expressed base transcript IDs from %s",
                len(expressed), counts_path)
    return expressed


def read_predictions(path: str) -> list:
    """Yield prediction dicts from the polyA_prediction.tsv file."""
    rows = []
    if path.endswith(".gz"):
        import gzip
        opener = lambda p: gzip.open(p, "rt")
    else:
        opener = open
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:len(TSV_COLUMNS)] != TSV_COLUMNS:
            raise SystemExit(
                "Unexpected polyA_prediction.tsv header. Expected %s, got %s"
                % (TSV_COLUMNS, header))
        for line in fh:
            if not line.strip():
                continue
            v = line.rstrip("\n").split("\t")
            rows.append({
                "chromosome": v[0],
                "transcript_id": v[1],
                "gene_id": v[2],
                "prediction": int(v[3]),
                "counts": int(v[4]),
                "flag": v[5],
            })
    return rows


def _prf(tp: int, n_pred: int, n_truth: int) -> tuple:
    precision = 100.0 * tp / n_pred if n_pred else 0.0
    recall = 100.0 * tp / n_truth if n_truth else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return precision, recall, f1


def compute_metrics(predictions: list, truth: dict, tolerance: int,
                    expressed: Optional[set] = None) -> dict:
    """Site-level evaluation. Each truth polyA site is an independent target.

    Precision: fraction of predictions that match some truth site.
    Recall:    fraction of truth sites matched by some prediction.
    Known/novel are tracked as subsets of truth sites.
    """
    if expressed is not None:
        truth = {tid: entries for tid, entries in truth.items()
                 if tid in expressed}
        predictions = [p for p in predictions if p["transcript_id"] in expressed]
        logger.info("After expression filter: %d truth transcripts, "
                    "%d predictions", len(truth), len(predictions))

    # Group predictions by transcript.
    preds_by_tid: dict = defaultdict(list)
    for p in predictions:
        preds_by_tid[p["transcript_id"]].append(p)

    # Count truth sites.
    n_sites_known = n_sites_novel = 0
    for entries in truth.values():
        for _, _, novel in entries:
            if novel:
                n_sites_novel += 1
            else:
                n_sites_known += 1
    n_sites = n_sites_known + n_sites_novel

    # Per-prediction: does it match any truth site?
    # The closest truth site determines known/novel category.
    tp_pred = fp_pred = unmatched = 0
    tp_pred_known = tp_pred_novel = 0
    fp_pred_known = fp_pred_novel = 0
    abs_distances: list = []
    for p in predictions:
        tid = p["transcript_id"]
        if tid not in truth:
            unmatched += 1
            continue
        entries = truth[tid]
        best_d: Optional[int] = None
        best_novel: bool = False
        for chr_id, pos, novel in entries:
            if chr_id != p["chromosome"]:
                continue
            d = abs(p["prediction"] - pos)
            if best_d is None or d < best_d:
                best_d = d
                best_novel = novel
        if best_d is None:
            unmatched += 1
        elif best_d <= tolerance:
            tp_pred += 1
            abs_distances.append(best_d)
            if best_novel:
                tp_pred_novel += 1
            else:
                tp_pred_known += 1
        else:
            fp_pred += 1
            if best_novel:
                fp_pred_novel += 1
            else:
                fp_pred_known += 1

    # Per-site: is the truth site matched by any prediction?
    recovered_known = 0
    recovered_novel = 0
    for tid, entries in truth.items():
        tid_preds = preds_by_tid.get(tid, [])
        for chr_id, pos, novel in entries:
            matched = any(
                p["chromosome"] == chr_id
                and abs(p["prediction"] - pos) <= tolerance
                for p in tid_preds)
            if matched:
                if novel:
                    recovered_novel += 1
                else:
                    recovered_known += 1
    recovered = recovered_known + recovered_novel

    n_predictions = len(predictions)
    precision = 100.0 * tp_pred / n_predictions if n_predictions else 0.0
    recall = 100.0 * recovered / n_sites if n_sites else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    mean_abs_distance = (sum(abs_distances) / len(abs_distances)
                         if abs_distances else 0.0)

    n_pred_known = tp_pred_known + fp_pred_known
    n_pred_novel = tp_pred_novel + fp_pred_novel
    precision_known = 100.0 * tp_pred_known / n_pred_known if n_pred_known else 0.0
    recall_known = 100.0 * recovered_known / n_sites_known if n_sites_known else 0.0
    f1_known = (2.0 * precision_known * recall_known / (precision_known + recall_known)
                if (precision_known + recall_known) > 0 else 0.0)
    precision_novel = 100.0 * tp_pred_novel / n_pred_novel if n_pred_novel else 0.0
    recall_novel = 100.0 * recovered_novel / n_sites_novel if n_sites_novel else 0.0
    f1_novel = (2.0 * precision_novel * recall_novel / (precision_novel + recall_novel)
                if (precision_novel + recall_novel) > 0 else 0.0)

    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "mean_abs_distance": mean_abs_distance,
        "n_predictions": n_predictions,
        "n_truth_sites": n_sites,
        "unmatched_predictions": unmatched,
        "tp": tp_pred,
        "fp": fp_pred,
        "precision_known": precision_known,
        "recall_known": recall_known,
        "f1_known": f1_known,
        "n_truth_known": n_sites_known,
        "precision_novel": precision_novel,
        "recall_novel": recall_novel,
        "f1_novel": f1_novel,
        "n_truth_novel": n_sites_novel,
    }


REPORT_KEYS = [
    "precision", "recall", "f1", "mean_abs_distance",
    "n_predictions", "n_truth_sites",
    "unmatched_predictions", "tp", "fp",
    "precision_known", "recall_known", "f1_known", "n_truth_known",
    "precision_novel", "recall_novel", "f1_novel", "n_truth_novel",
]


def write_report(metrics: dict, output_path: str) -> None:
    with open(output_path, "w") as fh:
        for k in REPORT_KEYS:
            v = metrics[k]
            if isinstance(v, float):
                fh.write("%s\t%.4f\n" % (k, v))
            else:
                fh.write("%s\t%d\n" % (k, v))


def main() -> int:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()
    if not os.path.exists(args.prediction):
        logger.error("Prediction file not found: %s", args.prediction)
        return 2
    if not os.path.exists(args.reference_gtf):
        logger.error("Reference GTF/DB not found: %s", args.reference_gtf)
        return 2

    db = get_or_build_db(args.reference_gtf)
    truth = load_truth(db)
    predictions = read_predictions(args.prediction)
    expressed = None
    if args.transcript_counts:
        if not os.path.exists(args.transcript_counts):
            logger.error("Transcript counts file not found: %s",
                         args.transcript_counts)
            return 2
        expressed = load_expressed_transcripts(args.transcript_counts)
    metrics = compute_metrics(predictions, truth, args.tolerance, expressed)
    write_report(metrics, args.output)

    logger.info("precision=%.4f recall=%.4f f1=%.4f n_pred=%d truth_sites=%d "
                "(known=%d novel=%d) recall_known=%.4f recall_novel=%.4f",
                metrics["precision"], metrics["recall"], metrics["f1"],
                metrics["n_predictions"], metrics["n_truth_sites"],
                metrics["n_truth_known"], metrics["n_truth_novel"],
                metrics["recall_known"], metrics["recall_novel"])
    logger.info("Report written to %s", args.output)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        print_exc()
        sys.exit(1)
