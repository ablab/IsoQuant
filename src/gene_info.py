############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import Enum, unique
from functools import partial
from collections import OrderedDict, defaultdict

from .serialization import *
from .common import (
    contains,
    equal_ranges,
    get_intron_strand,
    intervals_total_length,
    is_subprofile,
    overlaps,
    junctions_from_blocks,
    AtomicCounter
)

logger = logging.getLogger('IsoQuant')


@unique
class TranscriptModelType(Enum):
    known = 1
    novel_in_catalog = 2
    novel_not_in_catalog = 10


# simple class for storing all information needed for GFF
class TranscriptModel:
    def __init__(self, chr_id, strand, transcript_id, gene_id, exon_blocks, transcript_type):
        self.chr_id = chr_id
        self.strand = strand
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.exon_blocks = exon_blocks
        self.transcript_type = transcript_type
        self.additional_info = OrderedDict()
        self.intron_path = ()

    def get_start(self):
        return self.exon_blocks[0][0]

    def get_end(self):
        return self.exon_blocks[-1][1]

    def check_additional(self, attribute):
        return attribute in self.additional_info

    def add_additional_attribute(self, attribute, value):
        self.additional_info[attribute] = value

    def additional_attributes_str(self):
        return " ".join(['%s "%s";' % (k, v) for k,v in self.additional_info.items()])


# storage for feature profiles of all known isoforms of a gene or a set of overlapping genes
class FeatureProfiles:
    def __init__(self):
        self.profiles = {}
        self.profile_ranges = {}
        self.features = []

    def set_features(self, features):
        self.features = features

    def set_profiles(self, transcript_id, transcript_features, transcript_region, comaprator):
        self.profiles[transcript_id] = [-1] * len(self.features)
        current_profile = self.profiles[transcript_id]

        for pos, feature in enumerate(self.features):
            if not overlaps(feature, transcript_region):
                current_profile[pos] = -2

        pos = 0
        for feature in transcript_features:
            while pos < len(self.features) and not comaprator(feature, self.features[pos]):
                pos += 1
            while pos < len(self.features) and comaprator(feature, self.features[pos]):
                current_profile[pos] = 1
                pos += 1

        start_pos = 0
        while start_pos < len(current_profile) and current_profile[start_pos] == -2:
            start_pos += 1
        end_pos = len(current_profile) - 1
        while end_pos >= 0 and current_profile[end_pos] == -2:
            end_pos -= 1
        self.profile_ranges[transcript_id] = (start_pos, end_pos + 1)

    def print_debug(self):
        logger.debug(str(self.features))
        for t in self.profiles:
            logger.debug("Profile for isoform " + t)
            logger.debug(self.profiles[t])


# exon/intron info
class FeatureInfo:
    feature_id_counter = AtomicCounter()
    def __init__(self, chr_id, start, end, strand, type, gene_ids):
        self.id = FeatureInfo.feature_id_counter.increment()
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.strand = strand
        self.type = type
        self.gene_ids = gene_ids
        #self.id = "%s_%d_%d_%s" % (self.chr_id, self.start, self.end, self.strand)

    @staticmethod
    def header():
        return "#chr\tstart\tend\tstrand\tflags\tgene_ids"

    def to_str(self):
        return "%s\t%d\t%d\t%s\t%s\t%s" % (self.chr_id, self.start, self.end, self.strand, self.type, ",".join(self.gene_ids))


# All gene(s) information
class GeneInfo:
    EXTRA_BASES_FOR_SEQ = 20

    def __init__(self, gene_db_list, db, delta=0):
        if db is None:
            return
        # gffutils main structure
        self.db = db
        # list of genes in cluster
        self.gene_db_list = gene_db_list
        # gene region
        self.chr_id, self.start, self.end = self.get_gene_region()
        self.delta = delta

        # additional info for canonical splice site detection
        self.all_read_region_start = self.start
        self.all_read_region_end = self.end
        self.canonical_sites = {}
        self.gene_regions = {}
        self.reference_region = None

        # profiles for all known isoforoms
        self.intron_profiles = FeatureProfiles()
        self.exon_profiles = FeatureProfiles()
        self.split_exon_profiles = FeatureProfiles()

        self.all_isoforms_introns, self.all_isoforms_exons = self.set_introns_and_exons()
        self.split_exon_profiles.set_features(self.split_exons(self.exon_profiles.features))
        self.set_junction_profiles(self.all_isoforms_introns, self.all_isoforms_exons)

        self.isoform_strands = {}
        self.gene_strands = {}
        self.set_isoform_strands()
        self.gene_id_map = {}
        self.set_gene_ids()
        self.gene_attributes = {}
        self.set_gene_attributes()
        self.exon_property_map = self.set_feature_properties(self.all_isoforms_exons, self.exon_profiles)
        self.intron_property_map = self.set_feature_properties(self.all_isoforms_introns, self.intron_profiles)

    @classmethod
    def from_models(cls, transcript_model_storage, delta=0):
        gene_info = cls.__new__(cls)
        gene_info.db = None
        gene_info.gene_db_list = []
        if not transcript_model_storage:
            return cls([], None, delta)

        # gene region
        gene_info.chr_id = transcript_model_storage[0].chr_id
        gene_info.start = transcript_model_storage[0].get_start()
        gene_info.end = transcript_model_storage[0].get_end()
        gene_info.delta = delta
        gene_info.all_isoforms_exons = {}
        gene_info.all_isoforms_introns = {}
        gene_info.isoform_strands = {}
        gene_info.gene_id_map = {}
        gene_info.gene_attributes = {}
        introns = set()
        exons = set()

        for transcript_model in transcript_model_storage:
            gene_info.start = min(gene_info.start, transcript_model.get_start())
            gene_info.end = max(gene_info.end, transcript_model.get_end())
            t_id = transcript_model.transcript_id
            gene_info.all_isoforms_exons[t_id] = transcript_model.exon_blocks
            gene_info.all_isoforms_introns[t_id] = junctions_from_blocks(transcript_model.exon_blocks)
            exons.update(gene_info.all_isoforms_exons[t_id])
            introns.update(gene_info.all_isoforms_introns[t_id])
            gene_info.isoform_strands[transcript_model.transcript_id] = transcript_model.strand
            gene_info.gene_id_map[transcript_model.transcript_id] = transcript_model.gene_id

        # profiles for all known isoforoms
        gene_info.intron_profiles = FeatureProfiles()
        gene_info.exon_profiles = FeatureProfiles()
        gene_info.split_exon_profiles = FeatureProfiles()
        gene_info.ambiguous_isoforms = set()

        gene_info.intron_profiles.set_features(sorted(list(introns)))
        gene_info.exon_profiles.set_features(sorted(list(exons)))
        gene_info.split_exon_profiles.set_features(gene_info.split_exons(gene_info.exon_profiles.features))

        for transcript_model in transcript_model_storage:
            transcript_region = (transcript_model.get_start(), transcript_model.get_end())
            t_id = transcript_model.transcript_id
            gene_info.intron_profiles.set_profiles(t_id, gene_info.all_isoforms_introns[t_id], transcript_region, partial(equal_ranges, delta=0))
            gene_info.exon_profiles.set_profiles(t_id, gene_info.all_isoforms_exons[t_id], transcript_region, partial(equal_ranges, delta=0))
            gene_info.split_exon_profiles.set_profiles(t_id, gene_info.all_isoforms_exons[t_id], transcript_region, contains)

        gene_info.regions_for_bam_fetch = [(gene_info.start, gene_info.end)]
        gene_info.exon_property_map = None
        gene_info.intron_property_map = None

        # additional info for canonical splice site detection
        gene_info.all_read_region_start = gene_info.start
        gene_info.all_read_region_end = gene_info.end
        gene_info.canonical_sites = {}
        gene_info.gene_regions = {}
        gene_info.reference_region = None

        return gene_info

    @classmethod
    def from_model(cls, transcript_model, delta=0):
        gene_info = cls.__new__(cls)
        gene_info.db = None
        gene_info.gene_db_list = []
        # gene region
        gene_info.chr_id = transcript_model.chr_id
        gene_info.start = transcript_model.get_start()
        gene_info.end = transcript_model.get_end()
        gene_info.delta = delta

        # profiles for all known isoforoms
        gene_info.intron_profiles = FeatureProfiles()
        gene_info.exon_profiles = FeatureProfiles()
        gene_info.split_exon_profiles = FeatureProfiles()
        gene_info.ambiguous_isoforms = set()

        exons = transcript_model.exon_blocks
        introns = junctions_from_blocks(transcript_model.exon_blocks)
        t_id = transcript_model.transcript_id
        transcript_region = (transcript_model.get_start(), transcript_model.get_end())

        gene_info.all_isoforms_exons = {t_id: exons}
        gene_info.all_isoforms_introns = {t_id: introns}

        gene_info.exon_profiles.set_features(exons)
        gene_info.split_exon_profiles.set_features(exons)
        gene_info.intron_profiles.set_features(introns)

        gene_info.intron_profiles.set_profiles(t_id, introns, transcript_region, partial(equal_ranges, delta=0))
        gene_info.exon_profiles.set_profiles(t_id, exons, transcript_region, partial(equal_ranges, delta=0))
        gene_info.split_exon_profiles.set_profiles(t_id, exons, transcript_region, contains)

        gene_info.isoform_strands = {}
        gene_info.isoform_strands[transcript_model.transcript_id] = transcript_model.strand
        gene_info.gene_id_map = {}
        gene_info.gene_id_map[transcript_model.transcript_id] = transcript_model.gene_id
        gene_info.gene_attributes = {}

        gene_info.regions_for_bam_fetch = [(gene_info.start, gene_info.end)]
        gene_info.exon_property_map = None
        gene_info.intron_property_map = None

        # additional info for canonical splice site detection
        gene_info.all_read_region_start = gene_info.start
        gene_info.all_read_region_end = gene_info.end
        gene_info.canonical_sites = {}
        gene_info.gene_regions = {}
        gene_info.reference_region = None

        return gene_info

    @classmethod
    def from_region(cls, chr_id, start, end, delta=0, chr_record=None):
        gene_info = cls.__new__(cls)
        gene_info.db = None
        gene_info.gene_db_list = []
        # gene region
        gene_info.chr_id = chr_id
        gene_info.start = start
        gene_info.end = end
        gene_info.delta = delta

        # profiles for all known isoforms
        gene_info.intron_profiles = FeatureProfiles()
        gene_info.exon_profiles = FeatureProfiles()
        gene_info.split_exon_profiles = FeatureProfiles()
        gene_info.ambiguous_isoforms = set()

        gene_info.all_isoforms_exons = {}
        gene_info.all_isoforms_introns = {}
        gene_info.isoform_strands = {}
        gene_info.gene_id_map = {}
        gene_info.gene_attributes = {}
        gene_info.regions_for_bam_fetch = [(start, end)]
        gene_info.exon_property_map = None
        gene_info.intron_property_map = None

        # additional info for canonical splice site detection
        gene_info.all_read_region_start = gene_info.start
        gene_info.all_read_region_end = gene_info.end
        gene_info.canonical_sites = {}
        gene_info.gene_regions = {}
        gene_info.reference_region = None
        if chr_record:
            gene_info.reference_region = str(chr_record[gene_info.all_read_region_start - 1:gene_info.all_read_region_end + 1])

        return gene_info

    @classmethod
    def deserialize(cls, infile, genedb):
        gene_info = cls.__new__(cls)
        gene_info.db = genedb
        gene_info.delta = read_int(infile)
        gene_info.gene_db_list = []

        gene_count = read_int(infile)
        for i in range(gene_count):
            gene_id = read_string(infile)
            if gene_info.db:
                gene_info.gene_db_list.append(gene_info.db[gene_id])
        gene_info.chr_id = read_string(infile)
        gene_info.start = read_int(infile)
        gene_info.end = read_int(infile)

        gene_info.all_read_region_start = gene_info.start
        gene_info.all_read_region_end = gene_info.end

        # the rest is computed based on the database
        gene_info.reference_region = None
        gene_info.canonical_sites = {}
        gene_info.gene_regions = {}
        gene_info.intron_profiles = FeatureProfiles()
        gene_info.exon_profiles = FeatureProfiles()
        gene_info.split_exon_profiles = FeatureProfiles()
        gene_info.all_isoforms_introns = {}
        gene_info.all_isoforms_exons = {}
        if gene_info.gene_db_list:
            gene_info.all_isoforms_introns, gene_info.all_isoforms_exons = gene_info.set_introns_and_exons()
            gene_info.split_exon_profiles.set_features(gene_info.split_exons(gene_info.exon_profiles.features))
            gene_info.set_junction_profiles(gene_info.all_isoforms_introns, gene_info.all_isoforms_exons)

        gene_info.isoform_strands = {}
        gene_info.gene_strands = {}
        gene_info.set_isoform_strands()
        gene_info.gene_id_map = {}
        gene_info.set_gene_ids()
        gene_info.gene_attributes = {}
        gene_info.set_gene_attributes()
        gene_info.exon_property_map = gene_info.set_feature_properties(gene_info.all_isoforms_exons, gene_info.exon_profiles)
        gene_info.intron_property_map = gene_info.set_feature_properties(gene_info.all_isoforms_introns, gene_info.intron_profiles)
        return gene_info

    def serialize(self, outfile):
        write_int(self.delta, outfile)
        write_int(len(self.gene_db_list), outfile)
        for g in self.gene_db_list:
            write_string(g.id, outfile)
        write_string(self.chr_id, outfile)
        write_int(self.start, outfile)
        write_int(self.end, outfile)

    def empty(self):
        return not self.gene_db_list and not self.exon_profiles.features

    def print_debug(self):
        gene_names = []
        for g in self.gene_db_list:
            gene_names.append(g.id)
        logger.debug(">>> Processing gene set:" + ", ".join(gene_names))
        logger.debug("Split exon profiles:")
        self.split_exon_profiles.print_debug()
        logger.debug("Intron profiles:")
        self.intron_profiles.print_debug()
        logger.debug(">>> End of  gene informaiotn")

    # return region of overlapping gene set
    def get_gene_region(self):
        start = self.gene_db_list[0].start
        end = max(g.end for g in self.gene_db_list)
        chr_id = self.gene_db_list[0].seqid

        for gene_db in self.gene_db_list:
            if start > gene_db.start:
                start = gene_db.start
            if end < gene_db.end:
                end = gene_db.end

        return chr_id, start, end

    # return regions of every gene
    def get_gene_regions(self):
        if self.gene_regions:
            return self.gene_regions
        self.gene_regions = {}
        for gene_db in self.gene_db_list:
            self.gene_regions[gene_db.id] = (gene_db.start, gene_db.end)
        return self.gene_regions

    # set strands
    def set_isoform_strands(self):
        self.isoform_strands = {}
        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype=('transcript', 'mRNA')):
                self.isoform_strands[t.id] = t.strand
        self.gene_strands = {}
        for gene_db in self.gene_db_list:
            self.gene_strands[gene_db.id] = gene_db.strand

    # set isoform_id -> gene_id map
    def set_gene_ids(self):
        self.gene_id_map = {}
        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype=('transcript', 'mRNA')):
                self.gene_id_map[t.id] = gene_db.id

    def set_gene_attributes(self):
        self.gene_attributes = defaultdict(str)
        for gene_db in self.gene_db_list:
            for attr in gene_db.attributes.keys():
                if attr in ['gene_id', 'ID', 'level']:
                    continue
                self.gene_attributes[gene_db.id] += '%s "%s"; ' % (attr, gene_db.attributes[attr][0])

    # assigns an ordered list of all known exons and introns to self.exons and self.introns
    # returns 2 maps, isoform id -> intron / exon list
    def set_introns_and_exons(self):
        # dictionary: isoform id -> ordered list of intron coordinates
        all_isoforms_introns = {}
        # dictionary: isoform id -> ordered list of exon coordinates
        all_isoforms_exons = {}

        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype=('transcript', 'mRNA'), order_by='start'):
                all_isoforms_exons[t.id] = []
                for e in self.db.children(t, order_by='start'):
                    if e.featuretype == 'exon':
                        all_isoforms_exons[t.id].append((e.start, e.end))

                all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])

        if self.db and not all_isoforms_exons:
            logger.warning("Gene %s has no exons / transcripts, check your input annotation" % self.gene_db_list[0].id)

        introns = set()
        exons = set()
        for i in all_isoforms_introns.keys():
            introns.update(all_isoforms_introns[i])
        for i in all_isoforms_exons.keys():
            exons.update(all_isoforms_exons[i])

        self.intron_profiles.set_features(sorted(list(introns)))
        self.exon_profiles.set_features(sorted(list(exons)))

        return all_isoforms_introns, all_isoforms_exons

    # set feature properties for exon/intron counts
    def set_feature_properties(self, isoforms_to_feature_map, feature_profiles):
        similar_features = set()
        contained_features = set()
        # FIXME: change to interval tree instead of brute force
        for f1 in feature_profiles.features:
            for f2 in feature_profiles.features:
                if f1 == f2:
                    continue
                if equal_ranges(f1, f2, self.delta):
                    similar_features.add(f1)
                    similar_features.add(f2)
                if contains(f1, f2):
                    contained_features.add(f2)
                elif contains(f2, f1):
                    contained_features.add(f1)

        feature_to_isoform = defaultdict(list)
        for t in isoforms_to_feature_map.keys():
            isoform_features = isoforms_to_feature_map[t]
            if len(isoform_features) == 0:
                continue
            elif len(isoform_features) == 1:
                feature_to_isoform[isoform_features[0]].append((t, 'T'))
            else:
                feature_to_isoform[isoform_features[0]].append((t, 'T'))
                feature_to_isoform[isoform_features[-1]].append((t, 'T'))
                for e in isoform_features[1:-1]:
                    feature_to_isoform[e].append((t, ''))

        feature_properties = []
        for feature in feature_profiles.features:
            if all(x[1] == 'T' for x in feature_to_isoform[feature]):
                feature_type = "X"
            elif any(x[1] == 'T' for x in feature_to_isoform[feature]):
                feature_type = "T"
            else:
                feature_type = "I"

            if feature in similar_features:
                # similar features
                feature_type += "S"
            if feature in contained_features:
                # feature contained in anther one
                feature_type += "C"

            strands = set([self.isoform_strands[t[0]] for t in feature_to_isoform[feature]])
            strand_str = "".join(strands)
            gene_ids = set([self.gene_id_map[t[0]] for t in feature_to_isoform[feature]])

            if len(feature_to_isoform[feature]) == 1:
                # unique feature, appears only in one isoform
                feature_type += "U"
            elif len(gene_ids) > 1:
                # multiple genes
                feature_type += "M"

            feature_properties.append(FeatureInfo(self.chr_id, feature[0], feature[1], strand_str,
                                                  feature_type, list(gene_ids)))

        assert len(feature_properties) == len(feature_profiles.features)
        return feature_properties

    # split exons into non-overlapping covering blocks
    @staticmethod
    def split_exons(exons):
        exon_starts = sorted(map(lambda x: x[0], exons))
        exon_ends = sorted(map(lambda x: x[1], exons))

        current_state = 0
        starts_pos = 0
        ends_pos = 0
        last_border = -1
        exon_blocks = []

        while starts_pos < len(exon_starts):
            if exon_starts[starts_pos] <= exon_ends[ends_pos]:
                # do not consider the same exon star
                if starts_pos == 0 or exon_starts[starts_pos] > exon_starts[starts_pos - 1]:
                    cur_border = exon_starts[starts_pos]
                    if last_border != -1 and current_state > 0:
                        exon_blocks.append((last_border, cur_border - 1))
                    last_border = cur_border
                current_state += 1
                starts_pos += 1
            else:
                if last_border == -1 or current_state == 0:
                    print("Error, exon ends before the start")

                if ends_pos == 0 or exon_ends[ends_pos] > exon_ends[ends_pos - 1]:
                    cur_border = exon_ends[ends_pos]
                    exon_blocks.append((last_border, cur_border))
                    last_border = cur_border + 1
                current_state -= 1
                ends_pos += 1

        while ends_pos < len(exon_ends):
            if ends_pos == 0 or exon_ends[ends_pos] > exon_ends[ends_pos - 1]:
                cur_border = exon_ends[ends_pos]
                exon_blocks.append((last_border, cur_border))
                last_border = cur_border + 1
            current_state -= 1
            ends_pos += 1

        if current_state != 0:
            logger.critical("Unequal number of starts and ends")

        if exon_blocks != sorted(exon_blocks):
            logger.critical("Somehow block are unsorted")

        return exon_blocks

    # calculate junction profiles for known isoforms
    def set_junction_profiles(self, all_isoforms_introns, all_isoforms_exons):
        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype=('transcript', 'mRNA'), order_by='start'):
                transcript_region = self.transcript_region(t.id)
                # setting up intron profiles for current isoform
                self.intron_profiles.set_profiles(t.id, all_isoforms_introns[t.id], transcript_region,
                                                  partial(equal_ranges, delta=0))
                # setting up exon profiles for current isoform
                self.exon_profiles.set_profiles(t.id, all_isoforms_exons[t.id], transcript_region,
                                                partial(equal_ranges, delta=0))
                # setting up split exon profiles for current isoform
                self.split_exon_profiles.set_profiles(t.id, all_isoforms_exons[t.id], transcript_region, contains)

    def transcript_start(self, transcript_id):
        return self.all_isoforms_exons[transcript_id][0][0]

    def transcript_end(self, transcript_id):
        return self.all_isoforms_exons[transcript_id][-1][1]

    def transcript_region(self, transcript_id):
        return (self.transcript_start(transcript_id), self.transcript_end(transcript_id))

    def transcript_exon_count(self, transcript_id):
        return sum([1 if e == 1 else 0 for e in self.exon_profiles.profiles[transcript_id]])

    def total_transcript_length(self, transcript_id):
        return intervals_total_length(self.all_isoforms_exons[transcript_id])

    # closed interval
    def get_ref_seq(self, ref_start, ref_end):
        assert self.reference_region
        if ref_start < self.all_read_region_start or ref_end > self.all_read_region_end:
            logger.warning("Trying to extract ref seq outside of reference region: %d:%d" % (ref_start, ref_end))
        left_pos = ref_start - self.all_read_region_start
        right_pos = ref_end - self.all_read_region_start
        return self.reference_region[left_pos:right_pos+1]

    def set_reference_sequence(self, start, end, chr_record):
        self.all_read_region_start = start
        self.all_read_region_end = end
        self.reference_region = \
            str(chr_record[self.all_read_region_start - 1:self.all_read_region_end + 1])
        self.canonical_sites = {}


class StrandDetector:
    def __init__(self, chr_record):
        self.strand_dict = {}
        self.chr_record = chr_record

    def set_strand(self, intron, strand=None):
        if strand:
            self.strand_dict[intron] = strand
        elif self.chr_record:
            self.strand_dict[intron] = get_intron_strand(intron, self.chr_record)

    def get_strand(self, introns, has_polya=False, has_polyt=False):
        count_fwd = 0
        count_rev = 0
        for intron in introns:
            if not intron in self.strand_dict:
                strand = get_intron_strand(intron, self.chr_record)
                self.strand_dict[intron] = strand
            else:
                strand = self.strand_dict[intron]
            if strand == '+':
                count_fwd += 1
            elif strand == '-':
                count_rev += 1
        if count_fwd == count_rev:
            if has_polya and not has_polyt:
                return '+'
            elif has_polyt and not has_polya:
                return '-'
            return '.'
        return '+' if count_rev < count_fwd else '-'