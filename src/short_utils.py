import pysam
import gffutils

from .common import junctions_from_blocks

def get_introns(f, chromosome, start, end):
		samfile = pysam.AlignmentFile(f, "rb") 
		intr = samfile.find_introns(samfile.fetch(chromosome, start = start, stop = end))
		samfile.close()
		i_list = set()
		for i in intr.keys():
			if(type(i)!="int"): 
				i_list.add((i[0]+1,i[1]))
		return i_list

def ref_from_region(reference, region): 
	
	gene_list = list(reference.features_of_type('gene', limit = region, order_by=('seqid', 'start')))
	
	all_isoforms_introns = {}
	all_isoforms_exons = {}
	
	for gene_db in gene_list:
		for t in reference.children(gene_db, featuretype=('transcript', 'mRNA'), order_by='start'):
			all_isoforms_exons[t.id] = []
			for e in reference.children(t, order_by='start'):
				if e.featuretype == 'exon':
					all_isoforms_exons[t.id].append((e.start, e.end))

			all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])
	
	introns = set()
	for i in all_isoforms_introns.keys():
		introns.update(all_isoforms_introns[i])
		
	return introns
