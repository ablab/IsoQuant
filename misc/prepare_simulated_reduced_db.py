#!/usr/bin/env python3
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Split a simulation GTF into the reduced-db trio used by
``reduced_db_gffcompare.py`` for transcript-discovery assessment of
**alternative TSS/polyA**.

The simulation GTF is the full reference annotation (plain transcript IDs)
plus alternative-end variants whose ID is ``<base>_<coord>`` (same intron
chain, different terminal exon). Mapping alt-end variants onto the
reduced-db "novel" split lets the existing 3-terminal-delta assessment
measure how well alternative ends are recovered:

  <prefix>.expressed.gtf       all expressed transcripts          (-> full)
  <prefix>.expressed_kept.gtf  expressed plain-ID transcripts     (-> known)
  <prefix>.excluded.gtf        expressed alt-end variants         (-> novel)

"Expressed" = count > 0 in the simulated counts TSV. Alt-end IDs are matched
by ``_<digits>`` at the end (the simulation appends the alternative terminal
coordinate); plain GENCODE IDs (incl. ``_PAR_Y``) do not match.
"""
import argparse
import re
import sys
from traceback import print_exc

ALT_END = re.compile(r'_\d+$')


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sim_gtf", required=True, help="simulation GTF (full annotation + alt-end variants)")
    p.add_argument("--counts", required=True, help="simulated counts TSV (transcript_id<TAB>counts[...])")
    p.add_argument("--output_prefix", required=True, help="output prefix for the .expressed/.expressed_kept/.excluded GTFs")
    return p.parse_args()


def load_expressed(counts_path):
    expressed = set()
    for line in open(counts_path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        try:
            if float(f[1]) > 0:
                expressed.add(f[0])
        except ValueError:
            continue
    return expressed


def get_transcript_id(attrs):
    i = attrs.find('transcript_id "')
    if i < 0:
        return None
    return attrs[i + 15:].split('"', 1)[0]


def main():
    args = parse_args()
    expressed = load_expressed(args.counts)
    out_all = open(args.output_prefix + ".expressed.gtf", "w")
    out_kept = open(args.output_prefix + ".expressed_kept.gtf", "w")
    out_excl = open(args.output_prefix + ".excluded.gtf", "w")

    tx = {"all": set(), "kept": set(), "excl": set()}
    for line in open(args.sim_gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("transcript", "exon"):
            continue
        # Some upstream synthetic GTFs (the lrgasp human ones) carry a bogus
        # extra column 8 (literal "") before the attributes; the attributes
        # are always the last field. Read from there and write a normalized
        # 9-column line so gffcompare parses it (9-col input is unchanged).
        attrs = f[-1]
        tid = get_transcript_id(attrs)
        if tid is None or tid not in expressed:
            continue
        out_line = "\t".join(f[:8] + [attrs]) + "\n"
        out_all.write(out_line)
        is_tx = f[2] == "transcript"
        if is_tx:
            tx["all"].add(tid)
        if ALT_END.search(tid):
            out_excl.write(out_line)
            if is_tx:
                tx["excl"].add(tid)
        else:
            out_kept.write(out_line)
            if is_tx:
                tx["kept"].add(tid)

    for fh in (out_all, out_kept, out_excl):
        fh.close()
    print("Expressed transcripts written: all=%d  kept/plain(known)=%d  excluded/alt-end(novel)=%d"
          % (len(tx["all"]), len(tx["kept"]), len(tx["excl"])))


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception:
        print_exc()
        sys.exit(-1)
