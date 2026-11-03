#!/usr/bin/env python3
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Simulate full-length, polyA-free (FLNC-like) long reads from a transcript
expression profile with badread.

Each expressed transcript becomes a badread contig with ``depth=<count>``. The
fragment length is set far above any transcript length, so badread returns
every transcript whole (no truncation) and the number of reads per transcript
is proportional to its count. Reads are oriented 5'->3' like isoseq refine
FLNC reads, and GENCODE transcript sequences carry no polyA tail.

Outputs:
  <prefix>.fastq.gz          reads named <n>_<transcript_id>
  <prefix>.counts.tsv        transcript_id, simulated read count, CPM
  <prefix>.gene_counts.tsv   gene_id, simulated read count, CPM
  <prefix>.<kind>.gtf        truth GTFs (with --truth_gtf_prefix): the input
                             .expressed/.expressed_kept/.excluded GTFs reduced
                             to transcripts that received at least one read
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Tuple

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TRUTH_KINDS = ("expressed", "expressed_kept", "excluded")
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True,
                        help="expression profile TSV: transcript_id<TAB>count[<TAB>...]")
    parser.add_argument("--transcripts", required=True,
                        help="transcript FASTA (GENCODE headers: transcript_id|gene_id|...)")
    parser.add_argument("--reads", type=int, required=True, help="approximate number of reads to simulate")
    parser.add_argument("--output_prefix", required=True, help="output prefix")
    parser.add_argument("--threads", type=int, default=16, help="number of parallel badread runs [16]")
    parser.add_argument("--seed", type=int, default=1, help="base random seed, chunk i uses seed + i [1]")
    parser.add_argument("--identity", default="99.5,100,0.3",
                        help="badread --identity, mean,max,stdev [99.5,100,0.3]")
    parser.add_argument("--error_model", default="pacbio2021", help="badread error/qscore model [pacbio2021]")
    parser.add_argument("--badread", default="badread", help="badread executable [badread]")
    parser.add_argument("--truth_gtf_prefix",
                        help="prefix of .expressed/.expressed_kept/.excluded GTFs to reduce to simulated transcripts")
    return parser.parse_args()


def load_counts(path: str) -> Dict[str, float]:
    counts = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                count = float(fields[1])
            except ValueError:
                continue  # header
            if count > 0:
                counts[fields[0]] = count
    return counts


def load_transcripts(path: str, wanted: Dict[str, float]) -> Tuple[Dict[str, str], Dict[str, str]]:
    sequences: Dict[str, List[str]] = {}
    gene_of: Dict[str, str] = {}
    current = None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].strip().split("|")
                current = fields[0] if fields[0] in wanted else None
                if current is not None:
                    sequences[current] = []
                    gene_of[current] = fields[1] if len(fields) > 1 else fields[0]
            elif current is not None:
                sequences[current].append(line.strip())
    return {t: "".join(s) for t, s in sequences.items()}, gene_of


def write_depth_fasta(path: str, sequences: Dict[str, str], counts: Dict[str, float]) -> None:
    with open(path, "w") as out:
        for transcript_id in sorted(sequences):
            out.write(">%s depth=%g\n%s\n" % (transcript_id, counts[transcript_id], sequences[transcript_id]))


def run_badread(args, reference: str, chunk: int, quantity: int) -> str:
    fastq = "%s.chunk%02d.fastq" % (args.output_prefix, chunk)
    log = fastq + ".log"
    cmd = [args.badread, "simulate", "--reference", reference, "--quantity", str(quantity),
           "--length", "300000,0", "--identity", args.identity,
           "--error_model", args.error_model, "--qscore_model", args.error_model,
           "--seed", str(args.seed + chunk),
           "--start_adapter", "0,0", "--end_adapter", "0,0",
           "--junk_reads", "0", "--random_reads", "0", "--chimeras", "0", "--glitches", "0,0,0"]
    with open(fastq, "w") as out, open(log, "w") as err:
        subprocess.run(cmd, stdout=out, stderr=err, check=True)
    return fastq


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def collect_reads(chunk_files: List[str], output_fastq: str) -> Dict[str, int]:
    # badread header: @<uuid> <contig>,<+strand|-strand>,<start>-<end> length=...
    read_counts: Dict[str, int] = defaultdict(int)
    read_index = 0
    with gzip.open(output_fastq, "wt") as out:
        for chunk_file in chunk_files:
            with open(chunk_file) as handle:
                while True:
                    header = handle.readline()
                    if not header:
                        break
                    seq = handle.readline().rstrip("\n")
                    handle.readline()
                    qual = handle.readline().rstrip("\n")
                    transcript_id, strand = header.split()[1].split(",")[:2]
                    if strand == "-strand":
                        seq, qual = reverse_complement(seq), qual[::-1]
                    read_index += 1
                    read_counts[transcript_id] += 1
                    out.write("@%d_%s\n%s\n+\n%s\n" % (read_index, transcript_id, seq, qual))
    return read_counts


def write_counts(path: str, counts: Dict[str, int]) -> None:
    total = sum(counts.values())
    with open(path, "w") as out:
        for feature_id in sorted(counts):
            out.write("%s\t%.2f\t%.6f\n" % (feature_id, counts[feature_id], counts[feature_id] * 1e6 / total))


def reduce_truth_gtf(source: str, target: str, transcripts: Dict[str, int]) -> int:
    with open(source) as handle:
        lines = handle.readlines()
    kept_genes = set()
    for line in lines:
        match = TRANSCRIPT_ID_RE.search(line)
        if match and match.group(1) in transcripts:
            kept_genes.add(GENE_ID_RE.search(line).group(1))
    kept = 0
    with open(target, "w") as out:
        for line in lines:
            if line.startswith("#"):
                out.write(line)
                continue
            match = TRANSCRIPT_ID_RE.search(line)
            if match:
                if match.group(1) not in transcripts:
                    continue
                if line.split("\t")[2] == "transcript":
                    kept += 1
            else:
                gene_match = GENE_ID_RE.search(line)
                if not gene_match or gene_match.group(1) not in kept_genes:
                    continue
            out.write(line)
    return kept


def main() -> None:
    args = parse_args()
    counts = load_counts(args.counts)
    sequences, gene_of = load_transcripts(args.transcripts, counts)
    missing = len(counts) - len(sequences)
    if missing:
        print("Warning: %d expressed transcripts are absent from %s" % (missing, args.transcripts), file=sys.stderr)
    counts = {t: c for t, c in counts.items() if t in sequences}

    reference = args.output_prefix + ".depth.fasta"
    write_depth_fasta(reference, sequences, counts)
    # reads are whole transcripts drawn proportionally to counts, so the expected read length is
    # the count-weighted mean transcript length
    mean_length = sum(counts[t] * len(sequences[t]) for t in counts) / sum(counts.values())
    chunk_quantity = int(args.reads * mean_length / args.threads)
    print("Expressed transcripts: %d, count-weighted mean length %.0f, %d chunks x %d bp"
          % (len(counts), mean_length, args.threads, chunk_quantity), file=sys.stderr)

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        chunk_files = list(executor.map(lambda i: run_badread(args, reference, i, chunk_quantity),
                                        range(args.threads)))

    read_counts = collect_reads(chunk_files, args.output_prefix + ".fastq.gz")
    write_counts(args.output_prefix + ".counts.tsv", read_counts)
    gene_counts: Dict[str, int] = defaultdict(int)
    for transcript_id, count in read_counts.items():
        gene_counts[gene_of[transcript_id]] += count
    write_counts(args.output_prefix + ".gene_counts.tsv", gene_counts)
    print("Simulated %d reads from %d transcripts of %d genes"
          % (sum(read_counts.values()), len(read_counts), len(gene_counts)), file=sys.stderr)

    if args.truth_gtf_prefix:
        for kind in TRUTH_KINDS:
            kept = reduce_truth_gtf("%s.%s.gtf" % (args.truth_gtf_prefix, kind),
                                    "%s.%s.gtf" % (args.output_prefix, kind), read_counts)
            print("Truth %s: %d transcripts" % (kind, kept), file=sys.stderr)

    for chunk_file in chunk_files:
        os.remove(chunk_file)


if __name__ == "__main__":
    main()
