import argparse
import pysam
import gffutils

from src.illumina_exon_corrector import IlluminaExonCorrector
from src.alignment_info import AlignmentInfo
from src.common import junctions_from_blocks, overlaps
from src.short_utils import ref_from_region

def parse_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("--short", "-s", help="input sam/bam file containing alignment of the short reads", 
						type=str, dest = "short", required = True)
	parser.add_argument("--ont", "-i", help="input sam/bam file containing alignment of the long reads", 
						type=str, dest = "ont", required = True)
	parser.add_argument("--reference", "-r", help="input reference gtf file",
                        type=str, dest = "ref", required = True)
	parser.add_argument("--output", "-o", help="output folder location",
						type=str, dest = "out", required = True)
	return parser.parse_args()
	
def get_reference_db(ref_file, out_folder):
	if(ref_file[-2:] != "db"):
		gtf = ref_file
		args.output_exists = os.path.exists(out_folder)
		if not args.output_exists:
			os.makedirs(out_folder)
		ref = os.path.join(out_folder, "mouse_reference.db")
		gtf2db(gtf, ref, True)
	else:
		ref = ref_file
	return gffutils.FeatureDB(ref)


args = parse_args()
short_file = args.short
long_file = pysam.AlignmentFile(args.ont, "rb")
ref_db = get_reference_db(args.ref, args.out)
gene_list = list(ref_db.features_of_type('gene', order_by=('seqid', 'start')))
current_region = (0,0)
chromosome = 0
genes = []
before = 0
after = 0
for gene in gene_list:
	if overlaps((gene.start,gene.end), current_region) and gene.seqid == chromosome: # Chromosom mit abfragen?
		genes.append(gene)
		current_region = (current_region[0], max(current_region[1], gene.end))
	else:
		if genes :
			corrector = IlluminaExonCorrector(chromosome, current_region[0], current_region[1], short_file)
			reference_introns = ref_from_region(ref_db, (chromosome, current_region[0], current_region[1]))
			for alignment in long_file.fetch(chromosome, start = current_region[0], stop = current_region[1]):
				ai = AlignmentInfo(alignment)
				if not ai.read_exons:
					logger.warning("Read has no aligned exons")
					continue
				exons = ai.read_exons
				introns = set(junctions_from_blocks(exons))
				corrected_introns = corrector.correct_read(introns)
				before += len(introns.intersection(reference_introns))
				after += len(corrected_introns.intersection(reference_introns))
		genes = [gene]
		current_region = (gene.start, gene.end)
		chromosome = gene.seqid
if genes:
	corrector = IlluminaExonCorrector(chromosome, current_region[0], current_region[1], short_file)
	reference_introns = ref_from_region(ref_db, (chromosome, current_region[0], current_region[1]))
	for alignment in long_file.fetch(chromosome, start = current_region[0], stop = current_region[1]):
		ai = AlignmentInfo(alignment)
		if not ai.read_exons:
			logger.warning("Read has no aligned exons")
			continue
		exons = ai.read_exons
		introns = set(junctions_from_blocks(exons))
		corrected_introns = corrector.correct_read(introns)
		before += len(introns.intersection(reference_introns))
		after += len(corrected_introns.intersection(reference_introns))
		
print("Number of correct introns before correction:", before)
print("Number of correct introns after correction:", after)
