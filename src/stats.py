############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import logging
from collections import defaultdict

import pandas as pd

from src.common import proper_plural_form


logger = logging.getLogger('IsoQuant')


class EnumStats:
    def __init__(self):
        self.stats_dict = defaultdict(int)

    # element must be Enum
    def add(self, element):
        self.stats_dict[element] += 1

    def print_start(self, header_string = ""):
        if header_string:
            logger.info(header_string)
        for e in sorted(self.stats_dict.keys(), key=lambda x:x.name):
            logger.info("%s: %d" % (e.name, self.stats_dict[e]))


def transform_counts(path_to_csv, label):
    df = pd.read_csv(path_to_csv, sep='\t')
    df_features = df[:-3].copy()
    df_features.rename(columns={'count': label}, inplace=True)
    df_features.drop(columns='group_id', inplace=True)
    df_stats = df[-3:].copy()
    df_stats.rename(columns={'group_id': label}, inplace=True)
    df_stats.drop(columns='count', inplace=True)
    return pd.concat([df_features, df_stats], ignore_index=True)


def combine_counts(input_data, output):
    logger.info("Aggregating counts for " + proper_plural_form("sample", len(input_data.samples)))
    sample_0 = input_data.samples[0]
    gene_stats = transform_counts(sample_0.out_gene_counts_tsv, sample_0.label)
    transcript_stats = transform_counts(sample_0.out_transcript_counts_tsv, sample_0.label)

    for sample in input_data.samples[1:]:
        gene_stats = pd.merge(gene_stats, transform_counts(sample.out_gene_counts_tsv, sample.label),
                              on='feature_id', how='outer')
        transcript_stats = pd.merge(transcript_stats, transform_counts(sample.out_transcript_counts_tsv, sample.label),
                                    on='feature_id', how='outer')

    gene_stats.to_csv(os.path.join(output, 'combined_gene_counts.tsv'), sep='\t', index=False)
    transcript_stats.to_csv(os.path.join(output, 'combined_transcript_counts.tsv'), sep='\t', index=False)
    logger.info("Aggregation finished")
