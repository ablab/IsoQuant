import logging
import os
import subprocess
import json
import pysam
from src.common import get_path_to_program
from src.gtf2db import *
from src.input_data_storage import SampleData

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

        index = None if args.clean_start else find_stored_index(args)

        if index is None:
            index = index_reference(self.aligner, args)
            store_index(index, args)
        return index

    def map_reads(self, args):
        samples = []
        readable_names_dict = {}
        annotation_file = find_annotation(self.aligner, args)
        for sample in args.input_data.samples:
            bam_files = []
            for fastq_files in sample.file_list:
                for fastq_file in fastq_files:
                    bam_file = None if args.clean_start else find_stored_alignment(fastq_file, annotation_file, args)
                    if bam_file is None:
                        bam_file = align_fasta(self.aligner, fastq_file, annotation_file, args, sample.label, sample.aux_dir)
                        store_alignment(bam_file, fastq_file, annotation_file, args)
                    bam_files.append([bam_file])
                    if fastq_file in args.input_data.readable_names_dict:
                        readable_names_dict[bam_file] = args.input_data.readable_names_dict[fastq_file]

            samples.append(SampleData(bam_files, sample.label, sample.out_dir))
        args.input_data.samples = samples
        args.input_data.readable_names_dict = readable_names_dict
        args.input_data.input_type = "bam"
        return args.input_data


def get_aligner(aligner):
    path = get_path_to_program(aligner)
    if not path:
        logger.critical('{aligner} is not found! Make sure that {aligner} is in your PATH or use other alignment method'.format(aligner=aligner))
        exit(-1)
    return path


def find_stored_index(args):
    reference_filename = os.path.abspath(args.reference)

    with open(args.index_config_path, 'r') as f_in:
        converted_indexes = json.load(f_in)

    index_filename = converted_indexes.get(reference_filename, {}).get('index_filename')
    logger.debug('Searching for previously created index for {}'.format(reference_filename))
    if index_filename is None:
        return None
    index_mtime = converted_indexes.get(reference_filename, {}).get('index_mtime')
    reference_mtime = converted_indexes.get(reference_filename, {}).get('reference_mtime')
    kmer_size = converted_indexes.get(reference_filename, {}).get('kmer_size')
    if os.path.exists(reference_filename) and os.path.getmtime(reference_filename) == reference_mtime:
        if os.path.exists(index_filename) and os.path.getmtime(index_filename) == index_mtime:
            if KMER_SIZE[args.data_type] == kmer_size:
                logger.info('Index file found. Using {}'.format(index_filename))
                return index_filename
    return None


def store_index(index, args):
    reference_filename = os.path.abspath(args.reference)
    index = os.path.abspath(index)

    with open(args.index_config_path, 'r') as f_in:
        converted_indexes = json.load(f_in)
    converted_indexes[reference_filename] = {
        'index_filename': index,
        'reference_mtime': os.path.getmtime(reference_filename),
        'index_mtime': os.path.getmtime(index),
        'kmer_size': KMER_SIZE[args.data_type]
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
    elif aligner == "minimap2":
        if os.path.isfile(index_name):
            logger.debug('Reusing reference index ' + index_name)
            return index_name
        minimap2_path = get_aligner('minimap2')
        command = [minimap2_path, '-t', str(args.threads),
                   '-k', str(KMER_SIZE[args.data_type]),
                   '-w', '5', '-d', index_name, args.reference]
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
        if args.gtf('.db'):
            args.gtf = convert_db_to_gtf(args)
        return os.path.abspath(args.gtf)
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
                db2bed(bed_fname, args.genedb)
                store_bed(bed_fname, args)

        return os.path.abspath(bed_fname)


def align_fasta(aligner, fastq_file, annotation_file, args, label, out_dir):
    # TODO: fix paired end reads
    fastq_path = os.path.abspath(fastq_file)
    fname, ext = os.path.splitext(fastq_path.split('/')[-1])
    alignment_prefix = os.path.join(out_dir, label)

    prefix_name = fname
    if prefix_name.endswith(".fq") or prefix_name.endswith(".fastq") or prefix_name.endswith(".fa") or prefix_name.endswith(".fasta"):
        prefix_name, _ = os.path.splitext(prefix_name)
    hash_annotation = "_" + ("%x" % hash(annotation_file))[2:8] if annotation_file else ""
    hash_index = ("%x" % hash(args.index))[2:8]
    hash_fastq = ("%x" % hash(fastq_path))[2:8]

    alignment_bam_path = os.path.join(out_dir, label + '_' + prefix_name + '_%s_%s%s.bam' % (hash_fastq, hash_index, hash_annotation))
    logger.info("Aligning %s to the reference, alignments will be saved to %s" % (os.path.abspath(fastq_path),
                                                                                  os.path.abspath(alignment_bam_path)))
    alignment_sam_path = alignment_bam_path[:-4] + '.sam'

    log_fpath = os.path.join(args.output, "alignment.log")
    log_file = open(log_fpath, "a")
    # TODO: add star for barcoded reads
    if aligner == "starlong":
        star_path = get_aligner('STARlong')
        zcat_option = " --readFilesCommand zcat " if ext.endswith('gz') else ""

        # Simple
        # command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM
        #  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        annotation_opts = "" if not annotation_file else " --sjdbGTFfile " + annotation_file + " --sjdbOverhang 140 "
        command = '{exec_path} {zcat} --runThreadN {threads} --genomeDir {ref_index_name}  --readFilesIn {transcripts}  ' \
                  '--outSAMtype BAM Unsorted --outSAMattributes NH HI NM MD --outFilterMultimapScoreRange 1 ' \
                  '--outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen ' \
                  '-1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 ' \
                  ' --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 1000000 --seedPerWindowNmax 1000 '\
                  + annotation_opts + \
                  ' --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 ' \
                  '--outFileNamePrefix {alignment_out}'.format(exec_path=star_path,
                                                               zcat=zcat_option,
                                                               threads=str(args.threads),
                                                               ref_index_name=args.index,
                                                               transcripts=fastq_path,
                                                               alignment_out=alignment_prefix)
        logger.info("Running STAR (takes a while)")
        if subprocess.call(command.split(), stdout=log_file, stderr=log_file) != 0:
            logger.critical("STAR finished with errors! See " + log_fpath)
            exit(-1)
        logger.info("Sorting alignments")
        try:
            pysam.sort('-@', str(args.threads), '-o', alignment_bam_path, alignment_prefix + 'Aligned.out.bam')
        except pysam.SamtoolsError as err:
            logger.error(err.value)
            exit(-1)

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
        logger.info("Running minimap2 (takes a while)")
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
    except OSError as err:
        logger.error(err)
        exit(-1)

    return alignment_bam_path

