import pysam
import gffutils

from .common import junctions_from_blocks


def get_region_from_db(db, region): 
	
	gene_list = list(db.features_of_type('gene', limit = region, order_by=('seqid', 'start')))
	
	all_isoforms_introns = {}
	all_isoforms_exons = {}
	
	for gene_db in gene_list:
		for t in db.children(gene_db, featuretype=('transcript', 'mRNA'), order_by='start'):
			all_isoforms_exons[t.id] = []
			for e in db.children(t, order_by='start'):
				if e.featuretype == 'exon':
					all_isoforms_exons[t.id].append((e.start, e.end))

			all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])
	
	introns = set()
	for i in all_isoforms_introns.keys():
		introns.update(all_isoforms_introns[i])
		
	return introns
