############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Per-sample run summary.

IsoQuant computes plenty of QC numbers along the way, but most of them only ever
reached the log (alignment statistics, read assignment statistics, transcript model
statistics), while the barcode and UMI stages wrote small TSVs of their own. This
module gathers all of them into one object, which is then written next to the other
outputs as SAMPLE.summary_stats.json and rendered as SAMPLE.summary.html.

Nothing here recomputes anything expensive: in-memory counters are passed in by the
pipeline, everything else is read back from the stat files the run has just written,
so a resumed run produces the same report.
"""

import glob
import logging
import json
import os
from typing import Any, Dict, List, Optional, Set, Tuple

from isoquant_lib.assignment.isoform_assignment import ReadAssignmentType
from isoquant_lib.utils.file_utils import read_stats_tsv

logger = logging.getLogger('IsoQuant')

# A grouped counts file bigger than this is only scanned for its group count, not for
# per-group depth: the detailed pass is linear in the number of non-zero cells.
MAX_GROUPED_SCAN_BYTES = 2 * 1024 ** 3

# ... and the memory of that pass is linear in the number of groups instead, which the
# file size does not bound: a raw-barcode run has millions of them in a small file. Past
# this many the per-group maps are dropped and the groups are only counted.
MAX_GROUPED_SCAN_GROUPS = 1000000

# Keys of the UMI filtering stats file (UMIFilter in barcode_calling/umi_filtering.py):
# the reads that survived deduplication, and the reads that entered it.
UMI_SAVED_READS = "Total reads saved"
UMI_DEDUPLICATED_READS = "Assigned to any gene and barcoded"

# Feature whose per-group depth curve the report plots. The ranked list is one float
# per group, so it is only built for that one.
RANK_PLOT_FEATURE = "gene"


def enum_stats_to_dict(stats_dict: Dict[Any, int]) -> Dict[str, int]:
    """Convert an EnumStats counter ({enum member: count}) into a plain name -> count."""
    return {key.name: int(value) for key, value in sorted(stats_dict.items(), key=lambda kv: kv[0].name)}


# Headline assignment categories, in the order they are reported. Every assignment type
# ends up in exactly one of them, so the shares add up to the total; "other" holds the
# types that are none of the four (discarded and suspended reads).
ASSIGNMENT_BUCKETS = ("unique", "ambiguous", "inconsistent", "unassigned", "other")


def assignment_bucket(assignment_type: ReadAssignmentType) -> str:
    """Headline category of a single ReadAssignmentType.

    The enum's own predicates decide, so a type added to it later is classified rather
    than silently dropped. Inconsistent is checked before ambiguous because
    inconsistent_ambiguous and inconsistent_multigenic are both.
    """
    if assignment_type.is_unique():
        return "unique"
    if assignment_type.is_inconsistent():
        return "inconsistent"
    if assignment_type.is_unassigned():
        return "unassigned"
    if assignment_type.is_ambiguous():
        return "ambiguous"
    return "other"


def _rate(numerator: float, denominator: float) -> Optional[float]:
    if not denominator:
        return None
    return round(numerator / denominator, 4)


def _median(values: List[float]) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[middle])
    return (ordered[middle - 1] + ordered[middle]) / 2.0


class RunSummary:
    """Collects the QC numbers of a single experiment (sample)."""

    def __init__(self, sample_name: str, isoquant_version: str = "", command_line: str = "",
                 mode: str = ""):
        self.sample_name = sample_name
        self.isoquant_version = isoquant_version
        self.command_line = command_line
        self.mode = mode

        self.alignment: Dict[str, int] = {}
        # False when part of the input was not processed, which makes primary and
        # unaligned counts cover different sets of reads (see set_alignment_stats).
        self.alignment_covers_all_reads: bool = True
        self.assignment: Dict[str, int] = {}
        self.assignment_buckets: Dict[str, int] = {}
        self.transcript_models: Dict[str, int] = {}
        self.polya_reads: Optional[int] = None
        self.total_assignments: Optional[int] = None
        self.quantification: Dict[str, Dict[str, float]] = {}
        self.barcodes: Dict[str, int] = {}
        self.cell_barcodes: Dict[str, int] = {}
        self.umi_filtering: Dict[str, int] = {}
        self.umi_edit_distance: Optional[int] = None
        self.groups: Dict[str, object] = {}

    # ------------------------------------------------------------------ collectors

    def set_alignment_stats(self, stats_dict: Dict[Any, int],
                            covers_all_reads: bool = True) -> None:
        """stats_dict: EnumStats.stats_dict keyed by AlignmentType.

        covers_all_reads tells whether every reference carrying alignments was
        processed. It is not when --process_only_chr / --discard_chr are used, or when
        the annotation and the genome do not cover the same references: primary
        alignments are then counted for part of the input while the unaligned count
        comes from the whole file, and the totals are labelled as such.
        """
        self.alignment = enum_stats_to_dict(stats_dict)
        self.alignment_covers_all_reads = covers_all_reads

    def set_assignment_stats(self, stats_dict: Dict[ReadAssignmentType, int]) -> None:
        """stats_dict: EnumStats.stats_dict keyed by ReadAssignmentType.

        The headline categories are rolled up here, while the enum members are still
        at hand: their names alone do not say which category they belong to.
        """
        self.assignment = enum_stats_to_dict(stats_dict)
        rollup = {"total": 0}
        rollup.update({bucket: 0 for bucket in ASSIGNMENT_BUCKETS})
        for assignment_type, count in stats_dict.items():
            count = int(count)
            rollup["total"] += count
            rollup[assignment_bucket(assignment_type)] += count
        self.assignment_buckets = rollup

    def set_transcript_model_stats(self, stats_dict: Dict[Any, int]) -> None:
        """stats_dict: EnumStats.stats_dict keyed by TranscriptModelType."""
        self.transcript_models = enum_stats_to_dict(stats_dict)

    def set_polya_stats(self, total_assignments: int, polya_reads: int) -> None:
        self.total_assignments = total_assignments
        self.polya_reads = polya_reads

    def collect_output_files(self, sample: Any, grouping_strategy_names: Optional[List[str]] = None,
                             umi_edit_distance: Optional[int] = None) -> None:
        """Read back the stat files the run has written for this sample.

        umi_edit_distance is the filtering round whose survivors were counted; without
        it the lowest edit distance found on disk is described.
        """
        self._collect_counts_stats(sample)
        self._collect_barcode_stats(sample)
        self._collect_umi_stats(sample, umi_edit_distance)
        self._collect_group_stats(sample, grouping_strategy_names or [])

    def _collect_counts_stats(self, sample: Any) -> None:
        # dump_ungrouped() appends __ambiguous / __no_feature / __not_aligned to the
        # tail of the ungrouped counts file; everything above them is a real feature.
        for feature, counts_file in (("gene", getattr(sample, "out_gene_counts_tsv", None)),
                                     ("transcript", getattr(sample, "out_transcript_counts_tsv", None))):
            stats = self._read_counts_file(counts_file)
            if stats:
                self.quantification[feature] = stats

    @staticmethod
    def _read_counts_file(counts_prefix: Optional[str]) -> Dict[str, float]:
        if not counts_prefix:
            return {}
        counts_file = counts_prefix + "_counts.tsv"
        if not os.path.exists(counts_file):
            return {}
        assigned, special = 0.0, {}
        with open(counts_file) as f:
            next(f, None)  # header
            for line in f:
                values = line.rstrip("\n").split("\t")
                if len(values) != 2:
                    continue
                try:
                    count = float(values[1])
                except ValueError:
                    continue
                if values[0].startswith("__"):
                    special[values[0].lstrip("_")] = count
                else:
                    assigned += count
        if not assigned and not special:
            return {}
        # Only the counted reads are summed here: reads left out by the counting
        # strategy (inconsistent ones, typically) are in neither the features nor the
        # __ambiguous / __no_feature tail, so a ratio of the two would overstate the
        # assignment rate. The share of input reads is computed in to_dict() instead.
        return {"counted": assigned,
                "ambiguous": special.get("ambiguous", 0.0),
                "no_feature": special.get("no_feature", 0.0),
                "not_aligned": special.get("not_aligned", 0.0)}

    def _collect_barcode_stats(self, sample: Any) -> None:
        # One stats file per input file, written by the barcode calling stage.
        barcodes_tsv = getattr(sample, "barcodes_tsv", None)
        if barcodes_tsv:
            for stats_file in sorted(glob.glob(barcodes_tsv + "_*.tsv.stats")):
                read_stats_tsv(stats_file, self.barcodes)
        self.cell_barcodes = read_stats_tsv(getattr(sample, "out_cell_barcodes_stats", None))

    def _collect_umi_stats(self, sample: Any, edit_distance: Optional[int] = None) -> None:
        out_umi_filtered = getattr(sample, "out_umi_filtered", None)
        if not out_umi_filtered:
            return
        # <prefix>.UMI_filtered.ED<N>.stats.tsv, one per edit distance ever run here:
        # re-running the same directory in another mode leaves the previous one behind,
        # and its survivors are not the reads that were counted. Only the round the
        # caller names is reported; without one, the lowest edit distance is taken
        # (never the lexicographically last, which puts ED10 before ED2).
        if edit_distance is not None:
            candidates = [(edit_distance, "%s.ED%d.stats.tsv" % (out_umi_filtered, edit_distance))]
        else:
            candidates = sorted(self._umi_stat_files(out_umi_filtered))
        for found_distance, stats_file in candidates:
            stats = read_stats_tsv(stats_file)
            if not stats:
                continue
            self.umi_filtering = stats
            self.umi_edit_distance = found_distance
            return

    @staticmethod
    def _umi_stat_files(out_umi_filtered: str) -> List[Tuple[int, str]]:
        """(edit distance, path) for every UMI filtering stat file of this sample."""
        found = []
        for stats_file in glob.glob(out_umi_filtered + ".ED*.stats.tsv"):
            distance = os.path.basename(stats_file).split(".ED")[-1].split(".")[0]
            if distance.isdigit():
                found.append((int(distance), stats_file))
        return found

    def _collect_group_stats(self, sample: Any, grouping_strategy_names: List[str]) -> None:
        """Per-barcode/spot depth from the grouped counts written for each strategy."""
        for strategy in grouping_strategy_names:
            per_feature = {}
            for feature, prefix in (
                    ("gene", getattr(sample, "out_gene_grouped_counts_tsv", None)),
                    ("transcript", getattr(sample, "out_transcript_grouped_counts_tsv", None))):
                if not prefix:
                    continue
                linear_file = "%s_%s_counts.linear.tsv" % (prefix, strategy)
                stats = self._scan_linear_counts(linear_file,
                                                 keep_ranked=feature == RANK_PLOT_FEATURE)
                if stats:
                    per_feature[feature] = stats
            if per_feature:
                self.groups[strategy] = per_feature

    def _scan_linear_counts(self, linear_file: str, keep_ranked: bool = True) -> Dict[str, object]:
        """One sequential pass over a linear counts file: groups, reads per group and
        features per group. Returns {} when the file is absent.

        The detailed pass keeps two entries per group, so it gives up - on the file size
        up front and on the group count as it goes - and falls back to counting the
        groups, which needs one entry per group and no more.
        """
        if not os.path.exists(linear_file):
            return {}
        detailed = os.path.getsize(linear_file) <= MAX_GROUPED_SCAN_BYTES
        if not detailed:
            logger.info("%s is too large for per-group summary statistics, counting groups only"
                        % os.path.basename(linear_file))
        reads_per_group: Dict[str, float] = {}
        features_per_group: Dict[str, int] = {}
        groups: Set[str] = set()  # only filled once the detailed pass is given up
        total = 0.0
        with open(linear_file) as f:
            next(f, None)  # feature_id  group_id  count
            for line in f:
                values = line.rstrip("\n").split("\t")
                if len(values) != 3:
                    continue
                group_id = values[1]
                if not detailed:
                    groups.add(group_id)
                    continue
                try:
                    count = float(values[2])
                except ValueError:
                    continue
                total += count
                reads_per_group[group_id] = reads_per_group.get(group_id, 0.0) + count
                features_per_group[group_id] = features_per_group.get(group_id, 0) + 1
                if len(reads_per_group) > MAX_GROUPED_SCAN_GROUPS:
                    logger.info("%s holds more than %d groups, only counting them"
                                % (os.path.basename(linear_file), MAX_GROUPED_SCAN_GROUPS))
                    detailed = False
                    groups = set(reads_per_group)
                    reads_per_group, features_per_group, total = {}, {}, 0.0
        if not detailed:
            return {"groups": len(groups)} if groups else {}
        if not reads_per_group:
            return {}
        stats: Dict[str, object] = {"groups": len(reads_per_group), "reads": total}
        stats["median_reads_per_group"] = _median(list(reads_per_group.values()))
        stats["median_features_per_group"] = _median([float(v) for v in features_per_group.values()])
        if keep_ranked:
            # Depth per group, highest first: the barcode rank curve of the report.
            # Kept as a plain list, not the per-group dict, to stay small for datasets
            # with millions of barcodes.
            stats["ranked_reads"] = sorted(reads_per_group.values(), reverse=True)
        return stats

    # -------------------------------------------------------------------- derived

    @property
    def total_reads(self) -> Optional[int]:
        """Input reads: primary alignments plus the ones that did not align at all
        (secondary and supplementary alignments are extra records of the same reads).
        When part of the input was left out, primary alignments only cover the processed
        references, so this is the reads on processed references plus all unaligned
        reads, and the report labels it accordingly."""
        if not self.alignment:
            return None
        return self.alignment.get("primary", 0) + self.alignment.get("unaligned", 0)

    @property
    def mapping_rate(self) -> Optional[float]:
        return _rate(self.alignment.get("primary", 0), self.total_reads or 0)

    @property
    def barcode_rate(self) -> Optional[float]:
        return _rate(self.barcodes.get("Barcode detected", 0), self.barcodes.get("Total reads", 0))

    @property
    def umi_deduplicated_reads(self) -> Optional[int]:
        """Reads that entered UMI deduplication: assigned to a gene and barcoded.
        "Total assignments processed" also holds unbarcoded reads, which are never
        deduplicated and must not be counted as duplicates."""
        return self.umi_filtering.get(UMI_DEDUPLICATED_READS)

    @property
    def umi_duplication_rate(self) -> Optional[float]:
        """Share of the deduplicated reads that were removed as PCR/RT duplicates."""
        saved = self.umi_filtering.get(UMI_SAVED_READS)
        deduplicated = self.umi_deduplicated_reads
        if saved is None or not deduplicated:
            return None
        return _rate(deduplicated - saved, deduplicated)

    def assignment_rollup(self) -> Dict[str, int]:
        """Headline assignment categories on top of the per-type counts."""
        return dict(self.assignment_buckets)

    def group_rollup(self, strategy: str) -> Dict[str, object]:
        """Reads that carry a group id (barcode/spot) and were counted for a feature,
        as a share of all input reads - the single-cell "reads in cells" number."""
        per_feature = self.groups.get(strategy) or {}
        rollup = {}
        for feature, stats in per_feature.items():
            reads = stats.get("reads")
            rollup[feature] = {
                "groups": stats.get("groups"),
                "reads": reads,
                "share_of_reads": _rate(reads, self.total_reads or 0) if reads is not None else None,
                "median_reads_per_group": stats.get("median_reads_per_group"),
                "median_features_per_group": stats.get("median_features_per_group"),
            }
        return rollup

    # ----------------------------------------------------------------- serializing

    def to_dict(self) -> Dict[str, object]:
        summary: Dict[str, object] = {
            "sample": self.sample_name,
            "isoquant_version": self.isoquant_version,
            "mode": self.mode,
            "command_line": self.command_line,
        }
        if self.alignment:
            alignment = dict(self.alignment)
            alignment["total_reads"] = self.total_reads
            alignment["mapping_rate"] = self.mapping_rate
            alignment["covers_all_input_reads"] = self.alignment_covers_all_reads
            summary["alignment"] = alignment
        if self.assignment:
            assignment = {"by_type": dict(self.assignment)}
            assignment.update(self.assignment_rollup())
            if self.total_assignments is not None:
                assignment["total_assignments"] = self.total_assignments
                assignment["polya_reads"] = self.polya_reads
                assignment["polya_rate"] = _rate(self.polya_reads or 0, self.total_assignments)
            summary["assignment"] = assignment
        if self.quantification:
            quantification = {}
            for feature, stats in self.quantification.items():
                feature_stats = dict(stats)
                feature_stats["share_of_reads"] = _rate(stats.get("counted", 0), self.total_reads or 0)
                quantification[feature] = feature_stats
            summary["quantification"] = quantification
        if self.transcript_models:
            summary["transcript_models"] = dict(self.transcript_models)
        if self.barcodes:
            barcodes = {"by_stage": dict(self.barcodes), "barcode_rate": self.barcode_rate}
            if self.cell_barcodes:
                barcodes["cell_selection"] = dict(self.cell_barcodes)
            summary["barcodes"] = barcodes
        if self.umi_filtering:
            umi = {"by_stage": dict(self.umi_filtering), "edit_distance": self.umi_edit_distance}
            if self.umi_duplication_rate is not None:
                umi["molecules"] = self.umi_filtering[UMI_SAVED_READS]
                umi["duplication_rate"] = self.umi_duplication_rate
            summary["umi_filtering"] = umi
        if self.groups:
            summary["groups"] = {strategy: self.group_rollup(strategy) for strategy in self.groups}
        return summary

    def write_json(self, file_name: str) -> None:
        with open(file_name, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, sort_keys=False)
            f.write("\n")
