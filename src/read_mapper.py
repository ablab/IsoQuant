#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2023-2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
import shutil
import subprocess
import json
import pysam
import gzip

from .common import get_path_to_program
from .gtf2db import convert_db_to_gtf, db2bed
from .input_data_storage import SampleData

logger = logging.getLogger('IsoQuant')

PACBIO_CCS_DATA = 'pacbio_ccs'
NANOPORE_DATA = 'nanopore'
ASSEMBLY = 'assembly'

DATA_TYPE_ALIASES = {PACBIO_CCS_DATA: PACBIO_CCS_DATA, "pacbio": PACBIO_CCS_DATA, NANOPORE_DATA: NANOPORE_DATA,
                     "ont": NANOPORE_DATA, ASSEMBLY: ASSEMBLY, "transcripts": ASSEMBLY}
DATATYPE_TO_ALIGNER = {ASSEMBLY: 'minimap2', PACBIO_CCS_DATA: 'minimap2', NANOPORE_DATA: 'minimap2'}
                       # 'barcoded_se_reads' : 'star', 'barcoded_pe_reads' : 'star'}

SUPPORTED_ALIGNERS = ['starlong', 'minimap2']
SUPPORTED_STRANDEDNESS = ['forward', 'reverse', 'none']

KMER_SIZE = {ASSEMBLY: '15', PACBIO_CCS_DATA: '15', NANOPORE_DATA: '14'}
# CANONICAL_SITE_BONUS = {ASSEMBLY: '5', PACBIO_CCS_DATA: '5', NANOPORE_DATA: '9'}
MINIMAP_PRESET = {ASSEMBLY: 'splice:hq', PACBIO_CCS_DATA: 'splice:hq', NANOPORE_DATA: 'splice'}

class DataSetReadMapper:
    def __init__(self, args):
        self.args = args
        self.aligner = self.choose_aligner()
        self.index_fname = self.create_index(args)

    def choose_aligner(self):
        return self.args.aligner or DATATYPE_TO_ALIGNER[self.args.data_type]

    def create_index(self, args):
        if args.index and os.path.exists(args.index):
            return args.index

        index = None if args.clean_start else find_stored_index(args, self.aligner)

        if index is None:
            index = index_reference(self.aligner, args)
            store_index(index, args, self.aligner)
        return index

    def map_reads(self, args):
        samples = []
        annotation_file = find_annotation(self.aligner, args)
        for sample in args.input_data.samples:
            readable_names_dict = {}
            bam_files = []
            for fastq_files in sample.file_list:
                for fastq_file in fastq_files:
                    bam_file = None if args.clean_start else find_stored_alignment(fastq_file, annotation_file, args)
                    if bam_file is None:
                        bam_file = align_fasta(self.aligner, fastq_file, annotation_file, args, sample.prefix, sample.aux_dir)
                        store_alignment(bam_file, fastq_file, annotation_file, args)
                    bam_files.append([bam_file])
                    if fastq_file in sample.readable_names_dict:
                        readable_names_dict[bam_file] = sample.readable_names_dict[fastq_file]

            samples.append(SampleData(bam_files, sample.prefix, sample.out_dir, readable_names_dict, sample.illumina_bam))
        args.input_data.samples = samples
        args.input_data.input_type = "bam"
        return args.input_data


def get_aligner(aligner):
    path = get_path_to_program(aligner)
    if not path:
        logger.critical('{aligner} is not found! Make sure that {aligner} is in your PATH or use other alignment method'.format(aligner=aligner))
        exit(-1)
    return path


def find_stored_index(args, aligner):
    reference_filename = os.path.abspath(args.reference)

    with open(args.index_config_path, 'r') as f_in:
        converted_indexes = json.load(f_in)

    index_filename = converted_indexes.get(reference_filename, {}).get('index_filename')
    logger.debug('Searching for previously created index for {}'.format(reference_filename))
    if index_filename is None:
        return None
    index_mtime = converted_indexes.get(reference_filename, {}).get('index_mtime')
    reference_mtime = converted_indexes.get(reference_filename, {}).get('reference_mtime')
    used_aligner = converted_indexes.get(reference_filename, {}).get('aligner')
    if os.path.exists(reference_filename) and os.path.getmtime(reference_filename) == reference_mtime:
        if os.path.exists(index_filename) and os.path.getmtime(index_filename) == index_mtime:
            if aligner + str(KMER_SIZE[args.data_type]) == used_aligner:
                logger.info('Index file found. Using {}'.format(index_filename))
                return index_filename
    return None


def store_index(index, args, aligner):
    reference_filename = os.path.abspath(args.reference)
    index = os.path.abspath(index)

    with open(args.index_config_path, 'r') as f_in:
        converted_indexes = json.load(f_in)
    converted_indexes[reference_filename] = {
        'index_filename': index,
        'reference_mtime': os.path.getmtime(reference_filename),
        'index_mtime': os.path.getmtime(index),
        'aligner': aligner + str(KMER_SIZE[args.data_type])
    }
    with open(args.index_config_path, 'w') as f_out:
        json.dump(converted_indexes, f_out)
    logger.debug('New index saved to {}'.format(index))


def find_stored_bed(args):
    genedb_filename = os.path.abspath(args.genedb)

    with open(args.bed_config_path, 'r') as f_in:
        converted_beds = json.load(f_in)

    bed_filename = converted_beds.get(genedb_filename, {}).get('bed_filename')
    logger.debug('Searching for previously created BED for {}'.format(genedb_filename))
    if bed_filename is None:
        return None
    bed_mtime = converted_beds.get(genedb_filename, {}).get('bed_mtime')
    reference_mtime = converted_beds.get(genedb_filename, {}).get('reference_mtime')
    if os.path.exists(genedb_filename) and os.path.getmtime(genedb_filename) == reference_mtime:
        if os.path.exists(bed_filename) and os.path.getmtime(bed_filename) == bed_mtime:
                logger.info('BED file found. Using {}'.format(bed_filename))
                return bed_filename
    return None


def store_bed(bed, args):
    genedb_filename = os.path.abspath(args.genedb)
    bed = os.path.abspath(bed)

    with open(args.bed_config_path, 'r') as f_in:
        converted_beds = json.load(f_in)
    converted_beds[genedb_filename] = {
        'bed_filename': bed,
        'reference_mtime': os.path.getmtime(genedb_filename),
        'bed_mtime': os.path.getmtime(bed)
    }
    with open(args.bed_config_path, 'w') as f_out:
        json.dump(converted_beds, f_out)
    logger.debug('New BED saved to {}'.format(bed))


def find_stored_alignment(fastq_file, annotation, args):
    fastq = os.path.abspath(fastq_file)
    index = os.path.abspath(args.index)
    ann_path = os.path.abspath(annotation) if annotation else ""
    ann_str = "_" + ann_path if ann_path else ""

    key = "%s_aligned_to_%s%s" % (fastq, index, ann_str)
    with open(args.alignment_config_path, 'r') as f_in:
        aligned_fastq_files = json.load(f_in)

    logger.debug('Searching for previously created alignment for {}'.format(fastq))
    bam_fpath = aligned_fastq_files.get(key, {}).get('alignment_fpath')
    if bam_fpath is None:
        return None

    index_mtime = aligned_fastq_files.get(key, {}).get('index_mtime')
    if os.path.getmtime(index) != index_mtime:
        return None
    if ann_path:
        ann_mtime = aligned_fastq_files.get(key, {}).get('ann_mtime')
        if os.path.getmtime(ann_path) != ann_mtime:
            return None

    fastq_mtime = aligned_fastq_files.get(key, {}).get('fastq_mtime')
    bam_mtime = aligned_fastq_files.get(key, {}).get('bam_mtime')
    if os.path.exists(fastq) and os.path.getmtime(fastq) == fastq_mtime:
        if os.path.exists(bam_fpath) and os.path.getmtime(bam_fpath) == bam_mtime:
            logger.info('Bam alignment file found. Using {}'.format(bam_fpath))
            return bam_fpath
    return None


def store_alignment(bam_file, fastq_file, annotation, args):
    fastq = os.path.abspath(fastq_file)
    index = os.path.abspath(args.index)
    ann_path = os.path.abspath(annotation) if annotation else ""

    key = "%s_aligned_to_%s%s" % (fastq, index, "_" + ann_path if ann_path else "")
    bam_file = os.path.abspath(bam_file)

    with open(args.alignment_config_path, 'r') as f_in:
        aligned_fastq_files = json.load(f_in)
    aligned_fastq_files[key] = {
        'alignment_fpath': bam_file,
        'index_mtime': os.path.getmtime(index),
        'fastq_mtime': os.path.getmtime(fastq),
        'bam_mtime': os.path.getmtime(bam_file),
        'ann_mtime': os.path.getmtime(ann_path) if ann_path else ""
    }
    with open(args.alignment_config_path, 'w') as f_out:
        json.dump(aligned_fastq_files, f_out)
    logger.debug('New alignment saved to {}'.format(bam_file))


def index_reference(aligner, args):
    logger.info('Indexing reference')
    ref_name, _ = os.path.splitext(args.reference.split('/')[-1])

    command = ""
    index_name = os.path.join(os.path.abspath(args.output), "%s_k%s_idx" % (ref_name, KMER_SIZE[args.data_type]))

    if aligner == "starlong":
        if os.path.isdir(index_name) and os.path.exists(os.path.join(index_name, "genomeParameters.txt")):
            logger.debug('Reusing reference index ' + index_name)
            return index_name
        exec_path = get_aligner('STARlong')
        if not os.path.isdir(index_name):
            os.makedirs(index_name)
        command = [exec_path, '--runMode', 'genomeGenerate', '--runThreadN', str(args.threads),
                   '--genomeDir', index_name, '--genomeFastaFiles', args.reference]
        if args.indexing_options:
            command += args.indexing_options.split()

    elif aligner == "minimap2":
        if os.path.isfile(index_name):
            logger.debug('Reusing reference index ' + index_name)
            return index_name
        minimap2_path = get_aligner('minimap2')
        command = [minimap2_path, '-t', str(args.threads),
                   '-k', str(KMER_SIZE[args.data_type]),
                   '-w', '5']

        if args.indexing_options:
            command += args.indexing_options.split()
        command += ['-d', index_name, args.reference]

    else:
        logger.critical("Aligner " + aligner + " is not supported")
        exit(-1)

    log_fpath = os.path.join(args.output, "alignment.log")
    log_file = open(log_fpath, "w")
    if subprocess.call(command, stdout=log_file, stderr=log_file) != 0:
        logger.critical("Failed indexing reference! See " + log_fpath)
        exit(-1)
    return index_name


def find_annotation(aligner, args):
    if args.no_junc_bed:
        return None

    if aligner == "starlong":
        gene_annotation = args.genedb
        if gene_annotation.lower().endswith("db"):
            if args.original_annotation:
                gene_annotation = args.original_annotation
            else:
                gene_annotation = os.path.abspath(convert_db_to_gtf(args))

        if gene_annotation.lower().endswith("gz") or gene_annotation.lower().endswith("gzip"):
            genedb_file_name = os.path.basename(gene_annotation)
            genedb_name, _ = os.path.splitext(genedb_file_name)
            gunzipped_genedb = os.path.join(args.output, genedb_name)
            if not os.path.exists(gunzipped_genedb) or not args.resume:
                logger.info("Decompressing genome annotation to " + str(gunzipped_genedb))
                with open(gunzipped_genedb, "w") as outf:
                    shutil.copyfileobj(gzip.open(gene_annotation, "rt"), outf)
            gene_annotation = os.path.abspath(gunzipped_genedb)

        return os.path.abspath(gene_annotation)

    elif aligner == "minimap2":
        bed_fname = None
        if args.junc_bed_file:
            if not os.path.exists(args.junc_bed_file):
                logger.warning("Provided BED file %s does not exist, will generate from the annotation" %
                               args.junc_bed_file)
            else:
                bed_fname = args.junc_bed_file

        if bed_fname is None:
            bed_fname = find_stored_bed(args)
            if bed_fname is None:
                bed_fname = os.path.join(args.output, os.path.splitext(os.path.basename(args.genedb))[0] + ".bed")
                db2bed(args.genedb, bed_fname)
                store_bed(bed_fname, args)

        return os.path.abspath(bed_fname)


def align_fasta(aligner, fastq_file, annotation_file, args, label, out_dir):
    fastq_path = os.path.abspath(fastq_file)
    fname, ext = os.path.splitext(fastq_path.split('/')[-1])
    alignment_prefix = str(os.path.join(out_dir, label))

    prefix_name = fname
    if prefix_name.endswith(".fq") or prefix_name.endswith(".fastq") or prefix_name.endswith(".fa") or prefix_name.endswith(".fasta"):
        prefix_name, _ = os.path.splitext(prefix_name)
    hash_annotation = "_" + ("%x" % hash(annotation_file))[2:8] if annotation_file else ""
    hash_index = ("%x" % hash(args.index))[2:8]
    hash_fastq = ("%x" % hash(fastq_path))[2:8]

    alignment_bam_path = str(os.path.join(out_dir, label + '_' + prefix_name + '_%s_%s%s.bam' % (hash_fastq, hash_index, hash_annotation)))
    logger.info("Aligning %s to the reference, alignments will be saved to %s" % (os.path.abspath(fastq_path),
                                                                                  os.path.abspath(alignment_bam_path)))
    alignment_sam_path = alignment_bam_path[:-4] + '.sam'

    log_fpath = os.path.join(args.output, "alignment.log")
    log_file = open(log_fpath, "a")

    if aligner == "starlong":
        star_path = get_aligner('STARlong')
        zcat_option = ['--readFilesCommand', 'zcat'] if ext.endswith('gz') else []
        # command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM
        #  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        annotation_opts = [] if not annotation_file else ['--sjdbGTFfile', annotation_file, '--sjdbOverhang', '140']
        command = ([star_path, '--runThreadN', str(args.threads), '--genomeDir', args.index, '--readFilesIn',
                   fastq_path, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--seedPerReadNmax', '1000000',
                   '--outBAMsortingThreadN', str(args.threads), '--outSAMattributes', 'NH', 'HI', 'NM', 'MD'] +
                   annotation_opts + zcat_option + ['--outFileNamePrefix', alignment_prefix])
        if args.mapping_options:
            command += args.mapping_options.split()

        logger.info("Running STAR (takes a while)")
        if subprocess.call(command, stdout=log_file, stderr=log_file) != 0:
            logger.critical("STAR finished with errors! See " + log_fpath)
            exit(-1)
        shutil.move(alignment_prefix + "Aligned.sortedByCoord.out.bam", alignment_bam_path)

    elif aligner == "minimap2":
        additional_options = []
        minimap2_path = get_aligner('minimap2')
        if args.stranded == 'forward':
            additional_options.append('-uf')
        if annotation_file:
            args.junc_bed_file = annotation_file
            additional_options.append("--junc-bed")
            additional_options.append(annotation_file)

        command = [minimap2_path, args.index, fastq_path, '-a', '-x', MINIMAP_PRESET[args.data_type],
                   '--secondary=yes', '-Y', '--MD', '-t', str(args.threads)] + additional_options

        if args.mapping_options:
            command += args.mapping_options.split()

        minimap_version = "unknown"
        version_run = subprocess.run([minimap2_path, '--version'], capture_output=True)
        if version_run.returncode == 0:
            minimap_version = version_run.stdout.decode('UTF-8').strip()

        logger.info("Running minimap2 version %s (takes a while)" % minimap_version)
        if subprocess.call(command, stdout=open(alignment_sam_path, "w"), stderr=log_file) != 0:
            logger.critical("Minimap2 finished with errors! See " + log_fpath)
            exit(-1)
        logger.info("Sorting alignments")
        try:
            pysam.sort('-@', str(args.threads), '-o', alignment_bam_path, alignment_sam_path)
        except pysam.SamtoolsError as err:
            logger.error(err.value)
            exit(-1)

        os.remove(alignment_sam_path)
    else:
        logger.critical("Aligner " + aligner + " is not supported")
        exit(-1)
    logger.info("Indexing alignments")
    try:
        pysam.index(alignment_bam_path)
    except pysam.SamtoolsError as err:
        if "failed to create index" in err.value:
            logger.info(f"Samtools failed to generate default .bai index; with error: {err.value}")
            logger.info("Trying to build a CSI index instead")
            try:
                pysam.index('-@', str(args.threads), '-c', alignment_bam_path)
            except pysam.SamtoolsError as err:
                logger.error(f"Failed to create CSI index: {err.value}")
                exit(-1)
        else:
            logger.error(err.value)
            exit(-1) 
    except OSError as err:
        logger.error(err)
        exit(-1)

    return alignment_bam_path

