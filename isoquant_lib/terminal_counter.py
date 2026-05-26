############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Terminal position (polyA / TSS) prediction.

Per-transcript histogram of read end positions -> scipy peak detection ->
XGBoost peak filter (pre-trained model) -> classification against the
annotated transcript end. Predictions are written one row per accepted peak
(per group, when grouped) to a TSV alongside the standard counts files.

The shipped XGBoost classifiers live in ``isoquant_lib/data/``; they are
trained offline with ``scripts/train_polya_tss_model.py``.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import find_peaks, peak_prominences, peak_widths
from xgboost import XGBClassifier

from .isoform_assignment import (
    IsoformMatch,
    ReadAssignment,
    ReadAssignmentType,
)
from .long_read_counter import AbstractCounter
from .read_groups import AbstractReadGrouper

logger = logging.getLogger('IsoQuant')

_DATA_DIR = Path(__file__).parent / "data"
POLYA_MODEL_PATH = _DATA_DIR / "model_polya.json"
TSS_MODEL_PATH = _DATA_DIR / "model_tss.json"

# Peak finding + classification parameters. Match values used to train the
# shipped XGBoost models -- changing these without retraining will break the
# feature alignment between fit time and inference time.
PEAK_DISTANCE = 10
PEAK_REL_HEIGHT = 0.98
HISTOGRAM_PAD = 10
ANNOTATION_TOLERANCE = 10

_ACCEPTED_ASSIGNMENT_TYPES = frozenset({
    ReadAssignmentType.unique,
    ReadAssignmentType.unique_minor_difference,
    ReadAssignmentType.inconsistent,
    ReadAssignmentType.inconsistent_non_intronic,
})

# XGBoost feature columns (must match training)
FEATURE_COLUMNS = ['var', 'skew', 'peak_count', 'peak_width', 'entropy',
                   'mean_height', 'peak_heights', 'relative_height']

EMPTY_COLUMNS = ['chromosome', 'transcript_id', 'gene_id', 'prediction',
                 'counts', 'flag']
EMPTY_COLUMNS_GROUPED = EMPTY_COLUMNS + ['counts_byGroup', 'group_id']


class TerminalCounter(AbstractCounter):
    """Base class for polyA / TSS peak prediction counters.

    One counter instance handles a single chromosome (when constructed inside a
    worker) or the merged sample (when constructed at finalization). Each
    instance accumulates positional histograms per transcript across the reads
    it sees and emits its predictions in :meth:`dump`.
    """

    def __init__(self, args, output_prefix: str, model_path: Path,
                 string_pools=None, group_index: int = 0) -> None:
        # Skip AbstractCounter.__init__ -- we don't want counts_file_name's
        # suffix machinery; the prediction TSV path is already the full name.
        self.ignore_read_groups = string_pools is None
        self.output_prefix = output_prefix
        self.output_file = output_prefix
        self.output_counts_file_name = output_prefix
        self.output_tpm_file_name = None
        self.output_stats_file_name = None
        self.usable_file_name = None
        # Clear any stale per-chr file before we start appending in dump().
        open(self.output_file, "w").close()

        self.args = args
        self.string_pools = string_pools
        self.group_index = group_index

        # Model loading is deferred until dump() to keep XGBoost's OpenMP
        # runtime out of the parent process -- otherwise the inherited state
        # deadlocks fork() workers in ProcessPoolExecutor.
        self._model_path = model_path
        self._model = None

        # transcript_id -> {'chr', 'gene_id', 'data', 'annotated',
        #                   int_group_id -> list[int]}
        self.transcripts: dict = {}

    @property
    def model(self) -> XGBClassifier:
        if self._model is None:
            self._model = XGBClassifier()
            self._model.load_model(str(self._model_path))
        return self._model

    # -- AbstractCounter interface --------------------------------------------
    # Most of these are no-ops: the terminal counter only consumes finished
    # read assignments and emits its own TSV; it does not participate in the
    # confirmed-feature / unassigned / raw-feature flows.

    def add_read_info_raw(self, read_id, feature_ids, group_ids) -> None:
        return

    def add_confirmed_features(self, features) -> None:
        return

    def add_unassigned(self, read_assignment) -> None:
        return

    def add_unaligned(self, n_reads: int = 1) -> None:
        return

    def finalize(self, args=None) -> None:
        return

    def get_output_file_handler(self):
        return open(self.output_file, "a")

    # -- subclass hooks -------------------------------------------------------

    def _passes_filter(self, read_assignment: ReadAssignment) -> bool:
        raise NotImplementedError

    def _extract_position(self, read_assignment: ReadAssignment) -> Optional[int]:
        """Return the position of interest (genomic coordinate) or None."""
        raise NotImplementedError

    def _annotated_position(self, read_assignment: ReadAssignment,
                            transcript_id: str) -> int:
        raise NotImplementedError

    # -- accumulation ---------------------------------------------------------

    def add_read_info(self, read_assignment: Optional[ReadAssignment] = None) -> None:
        if read_assignment is None:
            return
        if read_assignment.assignment_type not in _ACCEPTED_ASSIGNMENT_TYPES:
            return
        if not self._passes_filter(read_assignment):
            return
        if not read_assignment.corrected_exons or not read_assignment.isoform_matches:
            return

        position = self._extract_position(read_assignment)
        if position is None:
            return

        isoform_match: IsoformMatch = read_assignment.isoform_matches[0]
        transcript_id = isoform_match.assigned_transcript
        if transcript_id is None:
            return

        entry = self.transcripts.get(transcript_id)
        if entry is None:
            entry = {
                'chr': read_assignment.chr_id,
                'gene_id': isoform_match.assigned_gene,
                'data': [],
                'annotated': self._annotated_position(read_assignment, transcript_id),
            }
            self.transcripts[transcript_id] = entry

        position = int(position)
        entry['data'].append(position)

        if not self.ignore_read_groups:
            group_id = self._read_group_id(read_assignment)
            entry.setdefault(group_id, []).append(position)

    def _read_group_id(self, read_assignment: ReadAssignment) -> int:
        if not read_assignment.read_group_ids:
            return 0
        return read_assignment.read_group_ids[self.group_index]

    # -- emission -------------------------------------------------------------

    def dump(self) -> None:
        if not self.transcripts:
            self._write_empty()
            return

        # Build base per-transcript dataframe.
        records = []
        for tid, entry in self.transcripts.items():
            data = entry['data']
            data_min = min(data)
            data_max = max(data)
            hist_counts = np.histogram(
                data,
                bins=(data_max - data_min) + 1,
                range=(data_min, data_max),
            )[0]
            padded = ([0] * HISTOGRAM_PAD + list(hist_counts) +
                      [0] * HISTOGRAM_PAD)
            records.append({
                'transcript_id': tid,
                'chromosome': entry['chr'],
                'gene_id': entry['gene_id'],
                'annotated': entry['annotated'],
                'start': data_min,
                'mode': int(stats.mode(data, keepdims=False).mode) - data_min,
                'mean': float(np.mean(data)) - data_min,
                'range': data_max - data_min + 1,
                'var': float(np.var(data)),
                'skew': (float(stats.skew(data))
                         if len(data) > 3 and np.var(data) >= 1e-12 else None),
                'entropy': float(stats.entropy(padded)),
                'mean_height': float(np.mean(padded)),
                'histogram': padded,
            })

        df = pd.DataFrame.from_records(records)

        # Peak detection per row. find_peaks returns positions relative to the
        # padded array; subtract HISTOGRAM_PAD to convert to start-relative.
        peaks_idx = df['histogram'].apply(
            lambda h: find_peaks(h, distance=PEAK_DISTANCE)[0])
        df['peak_count'] = peaks_idx.apply(len)
        df['peak_location'] = peaks_idx.apply(
            lambda p: [int(j - HISTOGRAM_PAD) for j in p])
        df['peak_prominence'] = df.apply(
            lambda r: peak_prominences(r['histogram'], peaks_idx.loc[r.name])[0],
            axis=1)
        widths = df.apply(
            lambda r: peak_widths(r['histogram'], peaks_idx.loc[r.name],
                                  rel_height=PEAK_REL_HEIGHT),
            axis=1)
        df['peak_width'] = widths.apply(lambda w: w[0])
        df['peak_heights'] = df.apply(
            lambda r: np.array([r['histogram'][int(p) + HISTOGRAM_PAD]
                                for p in r['peak_location']]),
            axis=1)
        df['peak_left'] = widths.apply(
            lambda w: [int(j - HISTOGRAM_PAD) for j in w[2]])
        df['peak_right'] = widths.apply(
            lambda w: [int(j - HISTOGRAM_PAD) for j in w[3]])

        # Transcripts with zero detected peaks: fall back to the mode.
        zero_peaks = df[df['peak_count'] == 0].copy()
        if not zero_peaks.empty:
            zero_peaks['prediction'] = zero_peaks['mode']
            zero_peaks['peak_location'] = zero_peaks['mode']
            zero_peaks['peak_heights'] = zero_peaks.apply(
                lambda r: r['histogram'][int(r['prediction']) + HISTOGRAM_PAD],
                axis=1)
            zero_peaks['peak_left'] = zero_peaks['mode']
            zero_peaks['peak_right'] = zero_peaks['mode']
            df = df[df['peak_count'] > 0].copy()

        # Transcripts with peaks: rank, explode, run XGBoost filter.
        if not df.empty:
            df = df.apply(self._rank_peaks, axis=1)
            df = df.explode(['peak_location', 'peak_prominence', 'peak_width',
                             'peak_heights', 'peak_left', 'peak_right',
                             'rank', 'relative_height']).reset_index(drop=True)
            df[FEATURE_COLUMNS] = df[FEATURE_COLUMNS].astype(float, errors='ignore')
            predicted = self.model.predict(df[FEATURE_COLUMNS].astype(float))
            df = df[predicted == 1].copy()
            df['prediction'] = df['peak_location']

        frames = [f for f in (zero_peaks, df) if not f.empty]
        if not frames:
            self._write_empty()
            return
        result = (pd.concat(frames, axis=0, ignore_index=True, sort=False)
                  if len(frames) > 1 else frames[0].reset_index(drop=True))

        # Convert peak position back to genomic coordinate and add flag.
        result['prediction'] = (result['prediction'].astype(int)
                                + result['start'].astype(int))
        result['flag'] = (
            (result['prediction'] - result['annotated']).abs()
                .gt(ANNOTATION_TOLERANCE).map({True: 'Novel', False: 'Known'}))
        result['counts'] = result.apply(self._counts_for_peak, axis=1)

        if self.ignore_read_groups:
            self._write_ungrouped(result)
        else:
            self._write_grouped(result)

    def _rank_peaks(self, row: pd.Series) -> pd.Series:
        heights = np.asarray(row['peak_heights'], dtype=float)
        if heights.size > 1:
            order = np.argsort(-heights)
            for key in ('peak_location', 'peak_prominence', 'peak_width',
                        'peak_heights', 'peak_left', 'peak_right'):
                row[key] = np.asarray(row[key])[order].tolist()
            row['rank'] = list(range(1, heights.size + 1))
            top = row['peak_heights'][0]
            row['relative_height'] = (
                [h / top if top else 0.0 for h in row['peak_heights']])
        else:
            row['rank'] = [0]
            row['relative_height'] = [1.0]
        return row

    def _counts_for_peak(self, row: pd.Series) -> int:
        hist = row['histogram']
        # peak_left / peak_right are start-relative; +HISTOGRAM_PAD to index
        # into the padded histogram. Clamp to its bounds in case of edge peaks.
        low = max(0, int(row['peak_left']) + HISTOGRAM_PAD)
        high = min(len(hist), int(row['peak_right']) + HISTOGRAM_PAD + 1)
        if high <= low:
            return 0
        return int(np.asarray(hist[low:high]).sum())

    def _counts_for_group(self, row: pd.Series, group_id: int) -> int:
        positions = self.transcripts[row['transcript_id']].get(group_id, [])
        if not positions:
            return 0
        # Count read ends falling within the genomic window of this peak.
        start = row['start']
        low = start + int(row['peak_left'])
        high = start + int(row['peak_right'])
        return int(sum(1 for p in positions if low <= p <= high))

    def _write_empty(self) -> None:
        cols = EMPTY_COLUMNS if self.ignore_read_groups else EMPTY_COLUMNS_GROUPED
        pd.DataFrame(columns=cols).to_csv(
            self.output_file, sep="\t", index=False, mode="w", header=True)

    def _write_ungrouped(self, result: pd.DataFrame) -> None:
        out = result[EMPTY_COLUMNS].copy()
        out.to_csv(self.output_file, sep="\t", index=False, mode="w", header=True)

    def _write_grouped(self, result: pd.DataFrame) -> None:
        # Collect all integer group IDs we have seen across the chromosome.
        seen_groups = sorted({gid for entry in self.transcripts.values()
                              for gid in entry.keys()
                              if isinstance(gid, int)})
        if not seen_groups:
            self._write_empty()
            return

        rows = []
        for _, peak in result.iterrows():
            for gid in seen_groups:
                count = self._counts_for_group(peak, gid)
                if count == 0:
                    continue
                rows.append({
                    'chromosome': peak['chromosome'],
                    'transcript_id': peak['transcript_id'],
                    'gene_id': peak['gene_id'],
                    'prediction': peak['prediction'],
                    'counts': peak['counts'],
                    'flag': peak['flag'],
                    'counts_byGroup': count,
                    'group_id': self._group_name(gid),
                })
        out = pd.DataFrame(rows, columns=EMPTY_COLUMNS_GROUPED)
        out.to_csv(self.output_file, sep="\t", index=False, mode="w", header=True)

    def _group_name(self, group_id: int) -> str:
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return pool.get_str(group_id)


class PolyACounter(TerminalCounter):
    """Predicts polyA cleavage sites per transcript."""

    def __init__(self, args, output_prefix: str,
                 string_pools=None, group_index: int = 0) -> None:
        super().__init__(args, output_prefix, POLYA_MODEL_PATH,
                         string_pools=string_pools, group_index=group_index)

    def _passes_filter(self, read_assignment: ReadAssignment) -> bool:
        return bool(read_assignment.polyA_found and read_assignment.polya_info)

    def _extract_position(self, read_assignment: ReadAssignment) -> Optional[int]:
        info = read_assignment.polya_info
        if read_assignment.strand == '+' and info.external_polya_pos != -1:
            return info.external_polya_pos
        if read_assignment.strand == '-' and info.external_polyt_pos != -1:
            return info.external_polyt_pos
        return None

    def _annotated_position(self, read_assignment: ReadAssignment,
                            transcript_id: str) -> int:
        exons = read_assignment.gene_info.all_isoforms_exons[transcript_id]
        if read_assignment.strand == '+':
            return exons[-1][1] + 1
        return exons[0][0] - 1


class TSSCounter(TerminalCounter):
    """Predicts transcription start sites per transcript."""

    def __init__(self, args, output_prefix: str,
                 string_pools=None, group_index: int = 0) -> None:
        super().__init__(args, output_prefix, TSS_MODEL_PATH,
                         string_pools=string_pools, group_index=group_index)

    def _passes_filter(self, read_assignment: ReadAssignment) -> bool:
        return read_assignment.strand in ('+', '-')

    def _extract_position(self, read_assignment: ReadAssignment) -> Optional[int]:
        exons = read_assignment.corrected_exons
        if not exons:
            return None
        if read_assignment.strand == '+':
            return exons[0][0]
        return exons[-1][1]

    def _annotated_position(self, read_assignment: ReadAssignment,
                            transcript_id: str) -> int:
        exons = read_assignment.gene_info.all_isoforms_exons[transcript_id]
        if read_assignment.strand == '+':
            return exons[0][0] - 1
        return exons[-1][1] + 1
