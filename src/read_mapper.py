import logging
import os
import subprocess

from Bio import SeqIO
from src.common import get_path_to_program
from src.gtf2db import convert_db_to_gtf
from src.input_data_storage import SampleData

logger = logging.getLogger('IsoQuant')

PACBIO_DATA = 'pacbio_raw'
PACBIO_CCS_DATA = 'pacbio_ccs'
NANOPORE_DATA = 'nanopore'
ASSEMBLY = 'assembly'

DATATYPE_TO_ALIGNER = {ASSEMBLY: 'minimap2', PACBIO_DATA: 'minimap2', PACBIO_CCS_DATA: 'minimap2',
                       NANOPORE_DATA: 'minimap2'}
                       # 'barcoded_se_reads' : 'star', 'barcoded_pe_reads' : 'star'}

SUPPORTED_ALIGNERS = ['starlong', 'minimap2']
SUPPORTED_STRANDEDNESS = ['forward', 'reverse', 'none']


class DataSetReadMapper:
    def __init__(self, args):
        self.args = args
        self.aligner = self.choose_aligner()
        self.index_fname = self.create_index(args)

    def choose_aligner(self,):
        if self.args.aligner is not None:
            return self.args.aligner
        else:
            return DATATYPE_TO_ALIGNER[self.args.data_type]

    def create_index(self, args):
        if args.index and os.path.exists(args.index):
            return args.index
        return index_reference(self.aligner, args)

    def map_reads(self, args):
        samples = []
        for sample in args.input_data.samples:
            bam_files = []
            for fastq_files in sample.file_list:
                bam_files.append([align_fasta(self.aligner, fastq_files, args)])
            samples.append(SampleData(bam_files, sample.label, sample.out_dir))
        args.input_data.samples = samples
        args.input_data.input_type = "bam"
        return args.input_data


def get_barcodes(contig_id, delim1, delim2):
    tokens = contig_id.split(delim1)
    if len(tokens) != 2:
        logger.warn("Wrong format " + contig_id)

    return tokens[0], tokens[1].strip().split(delim2) if len(tokens) > 1 else []


def convert_fasta_with_barcodes(contigs_file, args):
    fname, ext = os.path.splitext(contigs_file.split('/')[-1])
    short_id_contigs_fpath = (args.output, fname + "_short_ids.fasta")
    map_file_name = os.path.join(args.output, fname + "_map.txt")

    new_fasta = []
    with open(map_file_name, "w") as out_f:
        for record in SeqIO.parse(contigs_file, "fasta"):
            if len(record.seq) > args.max_len:
                continue
            name = record.id
            num = int(name.split('_')[1])
            if num % 2 == 1:
                continue
            record.id, bc = get_barcodes(name, args.delim, args.delim2)
            new_fasta.append(record)
            out_f.write(record.id + "_barcodeIDs_" + ",".join(bc) + "\n")

    SeqIO.write(new_fasta, short_id_contigs_fpath, "fasta")
    return short_id_contigs_fpath


def get_aligner(aligner):
    path = get_path_to_program(aligner)
    if not path:
        logger.critical('{aligner} is not found! Make sure that {aligner} is in your PATH or use other alignment method'.format(aligner=aligner))
        exit(-1)
    return path


def index_reference(aligner, args):
    logger.info('Indexing reference')
    ref_name, _ = os.path.splitext(args.reference.split('/')[-1])

    command = ""
    if aligner == "starlong":
        index_name = os.path.join(args.output, ref_name + "_idx")
        if os.path.isdir(index_name) and os.path.exists(os.path.join(index_name, "genomeParameters.txt")):
            logger.debug('Reusing reference index ' + index_name)
            return index_name
        # TODO: add annotation
        exec_path = get_aligner('STARlong')
        if not os.path.isdir(index_name):
            os.makedirs(index_name)
        command = [exec_path, '--runMode', 'genomeGenerate', '--runThreadN', str(args.threads),
                   '--genomeDir', index_name, '--genomeFastaFiles', args.reference]
    elif aligner == "minimap2":
        index_name = os.path.join(args.output, ref_name + ".idx")
        if os.path.isfile(index_name):
            logger.debug('Reusing reference index ' + index_name)
            return index_name
        minimap2_path = get_aligner('minimap2')
        k = 14 if args.data_type == NANOPORE_DATA else 15
        command = [minimap2_path, '-t', str(args.threads), '-k', str(k), '-w', '5', '-d', index_name, args.reference]
    else:
        logger.critical("Aligner " + aligner + " is not supported")
        exit(-1)

    log_fpath = os.path.join(args.output, "alignment.log")
    log_file = open(log_fpath, "w")
    if subprocess.call(command, stdout=log_file, stderr=log_file) != 0:
        logger.critical("Failed indexing reference! See " + log_fpath)
        exit(-1)
    return index_name


def align_fasta(aligner, fastq_paths, args):
    # TODO: fix paired end reads
    fastq_path = fastq_paths[0]
    logger.info("Aligning %s to the reference" % fastq_path)
    fname, ext = os.path.splitext(fastq_path.split('/')[-1])
    alignment_prefix = os.path.join(args.output, fname)
    alignment_sam_path = alignment_prefix + '.sam'
    alignment_bam_path = alignment_prefix + '.bam'

    # if 'barcoded' in args.data_type:
    #    fastq_path = convert_fasta_with_barcodes(fastq_path, args)

    log_fpath = os.path.join(args.output, "alignment.log")
    log_file = open(log_fpath, "a")
    # TODO: add star for barcoded reads
    if aligner == "starlong":
        star_path = get_aligner('STARlong')
        zcat_option = " --readFilesCommand zcat " if ext.endswith('gz') else ""

        # Simple
        # command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM
        #  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        gtf_fname = args.genedb
        if args.genedb.endswith('.db'):
            gtf_fname = convert_db_to_gtf(args)
        command = '{exec_path} {zcat} --runThreadN {threads} --genomeDir {ref_index_name}  --readFilesIn {transcripts}  ' \
                  '--outSAMtype BAM Unsorted --outSAMattributes NH HI NM MD --outFilterMultimapScoreRange 1 \
                   --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen ' \
                  '-1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 \
                   --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 1000000 --seedPerWindowNmax 1000 ' \
                  '--alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 ' \
                  '--sjdbGTFfile {annotation} --sjdbOverhang 140 \
                   --outFileNamePrefix {alignment_out}'.format(exec_path=star_path,
                                                               zcat=zcat_option,
                                                               threads=str(args.threads),
                                                               ref_index_name=args.index,
                                                               transcripts=fastq_path,
                                                               annotation=gtf_fname,
                                                               alignment_out=alignment_prefix)

        if subprocess.call(command.split(), stdout=log_file, stderr=log_file) != 0:
            logger.critical("STAR finished with errors! See " + log_fpath)
            exit(-1)
        subprocess.call(['samtools', 'sort', '-@', str(args.threads), '-o', alignment_bam_path,
                         alignment_prefix + 'Aligned.out.bam'], stderr=log_file)

    elif aligner == "minimap2":
        minimap2_path = get_aligner('minimap2')
        additional_options = []
        if args.stranded == 'forward':
            additional_options.append('-uf')
        # TODO: add junction bed
        command = [minimap2_path, args.index, fastq_path, '-a', '-x', 'splice', '--secondary=no', '-C', '5',
                   '-t', str(args.threads)] + additional_options
        if subprocess.call(command, stdout=open(alignment_sam_path, "w"), stderr=log_file) != 0:
            logger.critical("Minimap2 finished with errors! See " + log_fpath)
            exit(-1)
        subprocess.call(['samtools', 'sort', '-@', str(args.threads), '-o', alignment_bam_path, alignment_sam_path],
                        stderr=log_file)
    else:
        logger.critical("Aligner " + aligner + " is not supported")
        exit(-1)

    subprocess.call(['samtools', 'index', '-@', str(args.threads), alignment_bam_path])
    return alignment_bam_path

