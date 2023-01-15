#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import subprocess
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam
import gffutils
from enum import Enum, unique


@unique
class TranscriptType(Enum):
    known = 1
    novel = 2
    undefined = 0


class IsoQuantSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.find(".nic") != -1 or l.find(".nnic") != -1:
            return TranscriptType.novel
        elif l.find('transcript_id "SIRV') != -1 or l.find('transcript_id "ENS') != -1:
            return TranscriptType.known
        return TranscriptType.undefined


class StringTieSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.find("reference_id") != -1:
            return TranscriptType.known
        else:
            return TranscriptType.novel


class TranscriptIdSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.find('transcript_id "SIRV') != -1: # for SIRVs
            return TranscriptType.known
        elif l.find('transcript_id "ENS') != -1 and l.find("aligned_") == -1 and l.find("PB.") == -1: # for simulated data
            return TranscriptType.known
        else:
            return TranscriptType.novel


class FlamesSeparator:
    def __init__(self, _):
        pass

    def separate(self, l):
        if l.split('\t')[1] == "known":
            return TranscriptType.known
        else:
            return TranscriptType.novel


class CountTranscriptIdSeparator:
    def __init__(self, gtf_path):
        print("Reading counts")
        self.count_dict = defaultdict(float)
        for l in open(gtf_path + ".counts"):
            if l.startswith("#") or l.startswith("TXNAME"):
                continue
            t = l.strip().split()
            tid = t[0]
            self.count_dict[tid] = max(self.count_dict[tid], float(t[2]))

    def separate(self, l):
        tpos = l.find('transcript_id')
        if tpos == -1:
            return TranscriptType.undefined
        idpos = tpos + len('transcript_id') + 2
        endpos = l.find(";", idpos)
        if endpos == -1:
            print("Warning, unable to find ;")
            return TranscriptType.undefined
        tid = l[idpos:endpos-1]

        if tid not in self.count_dict or self.count_dict[tid] == 0:
            return TranscriptType.undefined
        elif tid.startswith('ENS') or tid.startswith('SIRV'):
            return TranscriptType.known
        else:
            return TranscriptType.novel


SEPARATE_FUNCTORS = {'isoquant':IsoQuantSeparator,
                     'stringtie':StringTieSeparator,
                     'flair':TranscriptIdSeparator,
                     'talon':TranscriptIdSeparator,
                     'bambu':CountTranscriptIdSeparator,
                     'flames':FlamesSeparator}


def split_gtf(ingtf_path, seaprator, out_full_path, out_known_path, out_novel_path):
    out_full = open(out_full_path, "w")
    out_known = open(out_known_path, "w")
    out_novel = open(out_novel_path, "w")
    for l in open(ingtf_path):
        if l.startswith("#"):
            continue
        ttype = seaprator.separate(l)
        if ttype == TranscriptType.undefined:
            continue
        out_full.write(l)
        if ttype == TranscriptType.novel:
            out_novel.write(l)
        elif ttype == TranscriptType.known:
            out_known.write(l)
    out_full.close()
    out_novel.close()
    out_known.close()


def run_gff_compare_noref(gtf_list, output):
    cmd_list = ["gffcompare", "-o", output] + gtf_list
    print(' '.join(cmd_list))
    result = subprocess.run(cmd_list)

    if result.returncode != 0:
        print("gffcompare failed ")
        return


def run_gff_compare(reference_gtf, compared_gtf, output, additional_option=""):
    cmd_list = ["gffcompare"]
    if additional_option:
        cmd_list.append(additional_option)
    cmd_list += ["-r", reference_gtf, "-o", output, compared_gtf]
    print(' '.join(cmd_list))
    result = subprocess.run(cmd_list)

    if result.returncode != 0:
        print("gffcompare failed for " + compared_gtf)
        return
