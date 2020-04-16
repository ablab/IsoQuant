import logging
import os

from Bio import SeqIO
from src.common import get_path_to_program

logger = logging.getLogger('IsoQuant')


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
    path = get_path_to_program(aligner.lower())
    if not path:
        logger.critical('{aligner} is not found! Make sure that {aligner} is in your PATH or use other alignment method'.format(aligner=aligner))
        exit(-1)
    return path


def index_reference(aligner, args):
    logger.debug('Indexing reference')
    ref_name, _ = os.path.splitext(args.reference.split('/')[-1])
    index_name = os.path.join(args.output, ref_name + ".idx")
    if os.path.exists(index_name):
        return index_name

    #TODO: fix
    if aligner == "starlong":
        pass
    elif aligner == "minimap2":
        minimap2_path = get_aligner('minimap2')
        command = '{exec_path} -t {threads} -k15 -w5 -d {ref_index_name} {ref_name}'.format(exec_path=minimap2_path,
                                                                                            threads=args.threads,
                                                                                            ref_index_name=index_name,
                                                                                            ref_name=args.reference)
    elif aligner == "gmap":
        pass

    if os.system(command) != 0:
        logger.critical("Failed indexing reference")
        exit(-1)
    return index_name

def align_fasta(aligner, fastq_paths, args):
    #TODO: fix paired end reads
    fastq_path = fastq_paths[0]
    fname, ext = os.path.splitext(fastq_path.split('/')[-1])
    alignment_prefix = os.path.join(args.output, fname)
    alignment_sam_path = alignment_prefix + '.sam'
    alignment_bam_path = alignment_prefix + '.bam'

    if 'barcoded' in args.data_type:
        fastq_path = convert_fasta_with_barcodes(fastq_path, args)

    # TODO: add star for barcoded reads
    if aligner == "starlong":
        #TODO: add indexing
        star_path = get_aligner('STARlong')
        zcat_option = " --readFilesCommand zcat " if ext.endswith('gz') else ""

        # Simple
        # command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM
        #  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        command = '{exec_path} {zcat} --runThreadN {threads} --genomeDir {ref_index_name}  --readFilesIn {transcripts}  ' \
                  '--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --outFilterMultimapScoreRange 1 \
                   --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen ' \
                  '-1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 \
                   --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 1000000 --seedPerWindowNmax 1000 ' \
                  '--alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 \
                   --outFileNamePrefix {alignment_out}'.format(exec_path=star_path,
                                                               zcat=zcat_option,
                                                               threads=str(args.threads),
                                                               ref_index_name=args.index,
                                                               transcripts=fastq_path,
                                                               alignment_out=alignment_prefix)

        exit_code = os.system(command)
        os.system('mv ' + alignment_prefix + 'Aligned.sortedByCoord.out.bam ' + alignment_bam_path)

        if exit_code != 0:
            logger.critical("STAR finished with errors")
            exit(-1)

    elif aligner == "minimap2":
        minimap2_path = get_aligner('minimap2')
        command = '{exec_path} {ref_index_name} {transcripts} -a -x splice --secondary=no -C5 -t {threads} ' \
                  '> {alignment_out} '.format(exec_path=minimap2_path,
                                              ref_index_name=args.index,
                                              transcripts=fastq_path,
                                              threads=str(args.threads),
                                              alignment_out=alignment_sam_path)
        exit_code = os.system(command)
        if exit_code != 0:
            logger.critical("Minimap2 finished with errors")
            exit(-1)
        os.system('samtools sort -o ' + alignment_bam_path + ' ' + alignment_sam_path)
    elif aligner == "gmap":
        #TODO: add indexing
        gmap_path = get_aligner('GMAP')
        command = gmap_path + ' -D ./  -d {ref_index_name} {transcripts} --format=samse -t {threads} -O ' \
                              '> {alignment_out} '.format(ref_index_name=args.index,
                                                          transcripts=fastq_path,
                                                          threads=str(args.threads),
                                                          alignment_out=alignment_sam_path)
        exit_code = os.system(command)
        if exit_code != 0:
            logger.critical("GMAP finished with errors")
            exit(-1)
        os.system('samtools sort -o ' + alignment_bam_path + ' ' + alignment_sam_path)
    else:
        logger.critical("Aligner " + aligner + " is not supported")
        exit(-1)

    os.system('samtools index ' + alignment_bam_path)
    return alignment_bam_path

