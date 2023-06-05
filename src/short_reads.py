import pysam
import sys
import argparse
import gffutils
import time
import networkx as nx
from networkx.algorithms import bipartite
from enum import Enum, unique

from common import *
from gtf2db import *

ABSENT_PAIR = (-1,-1)

@unique
class IntronType(Enum):
	equal = 1
	overlap = 4
	only_short = 2
	only_isoquant = 3
	short_contains_iso = 6
	iso_contains_short = 5
	
@unique
class IntronReference(Enum):
	short_right = 1
	iso_right = 2
	both = 3
	none = 4

# just doing everything here now to avoid messing with anything
def parse_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("--short", "-s", help="input sam/bam file containing alignment of the short reads", 
						type=str, dest = "short", required = True)
	parser.add_argument("--isoquant", "-i", help="input isoquant annotation", 
						type=str, dest = "iso", required = True)
	parser.add_argument("--reference", "-r", help="input reference gtf file",
                        type=str, dest = "ref", default = "")
	parser.add_argument("--output", "-o", help="output folder location",
						type=str, dest = "out", required = True)
	return parser.parse_args()

def get_introns(f):
	samfile = pysam.AlignmentFile(f, "rb") 
	introns = samfile.find_introns(read for read in samfile.fetch())
	samfile.close()
	i_list = set()
	for i in introns.keys():
		#i_list.append((i[0]+1,i[1]+1))
		i_list.add(i)
	return i_list
	

def introns_from_gene(db, gene, f):
	all_isoforms_introns = {}
	all_isoforms_exons = {}
	s = gene.start
	end = gene.end  # probably not just in the gene right? maybe also slightly before and slightly after?
	
	for t in db.children(gene, featuretype=('transcript', 'mRNA'), order_by='start'):
		all_isoforms_exons[t.id] = []
		for e in db.children(t, order_by='start'):
			if e.featuretype == 'exon':
				all_isoforms_exons[t.id].append((e.start, e.end))
				
		all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])
		
	introns = []
	for i in all_isoforms_introns.keys():
		introns.extend(all_isoforms_introns[i])
	
	samfile = pysam.AlignmentFile(f, "rb") 
	intr = samfile.find_introns(samfile.fetch(start = s, stop = end))
	samfile.close()
	i_list = []
	for i in intr.keys():
		#i_list.append((i[0]+1,i[1]+1))
		if(type(i)!="int"): 
			i_list.append((i[0]+1,i[1]))
		
	return introns, i_list
	
def introns_from_region(db, gene_list, current_region, f): # Kann ich das aus der reference auch selektiv rausbekommen? Vielleicht kann ich mir ne Liste ausgeben lassen innerhalb bestimmter Grenzen
	
	chromosome = gene_list[0].seqid
	
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
		
	samfile = pysam.AlignmentFile(f, "rb") 
	intr = samfile.find_introns(samfile.fetch(chromosome, start = current_region[0], stop = current_region[1]))
	samfile.close()
	i_list = set()
	for i in intr.keys():
		if(type(i)!="int"): 
			i_list.add((i[0]+1,i[1]))
	return introns, i_list, intr
	
def ref_from_region(reference, region): 
	
	gene_list = list(reference.features_of_type('gene', limit = region, order_by=('seqid', 'start')))
	
	all_isoforms_introns = {}
	all_isoforms_exons = {}
	
	for gene_db in gene_list:
		if gene_db.seqid == chromosome:
			for t in reference.children(gene_db, featuretype=('transcript', 'mRNA'), order_by='start'):
				all_isoforms_exons[t.id] = []
				for e in reference.children(t, order_by='start'):
					if e.featuretype == 'exon':
						all_isoforms_exons[t.id].append((e.start, e.end))

				all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])
	
	introns = set()
	for i in all_isoforms_introns.keys():
		introns.update(all_isoforms_introns[i])
		#print(i)
		#print(all_isoforms_introns[i])
		
	return introns

# from gene_info.py
def introns_from_db(db):
	gene_list = list(db.all_features(featuretype="gene"))
	#print(gene_list)
	all_isoforms_introns = {}
	all_isoforms_exons = {}
	for gene_db in gene_list:
		for t in db.children(gene_db, featuretype=('transcript', 'mRNA'), order_by='start'):
			all_isoforms_exons[t.id] = []
			for e in db.children(t, order_by='start'):
				if e.featuretype == 'exon':
					all_isoforms_exons[t.id].append((e.start, e.end))

			all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])

	if db and not all_isoforms_exons:
		logger.warning("Gene %s has no exons / transcripts, check your input annotation" % gene_list[0].id)
	
	introns = set()
	for i in all_isoforms_introns.keys():
		introns.update(all_isoforms_introns[i])
	
	return introns

def compare_introns(short, iso):
	classification = {}
	found_short = [False]*len(short)
	found_iso = [False]*len(iso)
	#print(len(short))
	#print(len(iso))
	for i in range(len(short)): 
		for j in range(len(iso)):
			#print(short[i][1])
			#print(iso[j][0])
			if short[i][1] < iso[j][0] or short[i][0] > iso[j][1]:
				#print("done")
				continue
			elif short[i] == iso[j]:
				classification[(short[i],iso[j])] = IntronType.equal
				found_short[i] = True
				found_iso[j] = True
				#print("equal")
				continue
			elif (short[i][0] <= iso[j][0] and short[i][1] <= iso[j][1]) or (short[i][0] >= iso[j][0] and short[i][1] >= iso[j][1]):
				classification[(short[i],iso[j])] = IntronType.overlap
				found_short[i] = True
				found_iso[j] = True
				#print("overlap")
				continue
			elif short[i][0] >= iso[j][0] and short[i][1] <= iso[j][1]:
				classification[(short[i],iso[j])] = IntronType.iso_contains_short
				found_short[i] = True
				found_iso[j] = True
				#print("short contained in iso")
				continue
			elif short[i][0] <= iso[j][0] and short[i][1] >= iso[j][1]:
				classification[(short[i],iso[j])] = IntronType.short_contains_iso
				found_short[i] = True
				found_iso[j] = True
				#print("iso contained in short")
				continue
		if not found_short[i]:
			classification[(short[i],ABSENT_PAIR)] = IntronType.only_short
			#print("not found in iso")
	for j in range(len(iso)):
		if not found_iso[j]:
			classification[(ABSENT_PAIR,iso[j])] = IntronType.only_isoquant
			#print("not found in short")
	return classification
	
def classify_introns(iso, short):
	if (short[0] <= iso[0] and short[1] <= iso[1]) or (short[0] >= iso[0] and short[1] >= iso[1]):
		return IntronType.overlap
	elif short[0] >= iso[0] and short[1] <= iso[1]:
		return IntronType.iso_contains_short
	elif short[0] <= iso[0] and short[1] >= iso[1]:
		return IntronType.short_contains_iso
	return IntronType.only_isoquant
	
def classify_reference(iso, short, ref, chromosome):
	#if iso in ref and any(sh in ref for sh in short):
	if iso in ref and short in ref:
		print("both:", chromosome, ":" , iso, ",", short)
		return IntronReference.both
	elif iso in ref:
		print("iso:", chromosome, ":" , iso, ",", short)
		return IntronReference.iso_right
	#elif any(sh in ref for sh in short):
	elif short in ref:
		print("short:", chromosome, ":" , iso, ",", short)
		return IntronReference.short_right
	else:
		print("none:", chromosome, ":" , iso, ",", short)
		return IntronReference.none
	
def compare_alg(short, iso, counts, ref, chromosome):
	equal = iso.intersection(short)
	short = short.difference(equal)
	iso = iso.difference(equal)
	G = nx.Graph()
	G.add_nodes_from(iso, bipartite = 0)
	G.add_nodes_from(short, bipartite = 1)
	score = 1000000000000
	sh = (0,0)
	#score = [1000000000000,1000000000000,1000000000000]
	#sh = [(0,0),(0,0),(0,0)]
	chosen = []
	classification = {}
	cl_ref = {}
	scounts = []
	for i in iso:
		for s in short: 
			x = abs(i[0] - s[0]) + abs(i[1] - s[1])
			if (overlaps(i, s) and abs(i[0] - s[0]) <= 12 and abs(i[1] - s[1]) <= 12):
				# if x < score[0]:
					# score[0] = x
					# sh[0] = s
				# elif x < score[1]:
					# score[1] = x
					# sh[1] = s
				# elif x < score[2]:
					# score[2] = x
					# sh[2] = s
				if x < score:
					score = x
					sh = s
				G.add_edge(i, s, weight = abs(i[0] - s[0]) + abs(i[1] - s[1]))
		#classification[(i,sh[0])] = classify_introns(i,sh[0])
		#classification[(i,sh[1])] = classify_introns(i,sh[1])
		#classification[(i,sh[2])] = classify_introns(i,sh[2])
		#classification[(i,sh)] = classify_introns(i,sh)
		classify_reference(i, sh, ref, chromosome)
		# for y in sh:
			# if(not y == (0,0)):
				# chosen.append(y)
				# scounts.append(counts[(y[0]-1,y[1])])
		if(not sh == (0,0)):
			chosen.append(sh)
			scounts.append(counts[(sh[0]-1,sh[1])])
		#score = [1000000000000,1000000000000,1000000000000]
		#sh = [(0,0),(0,0),(0,0)]
		score = 1000000000000
		sh = (0,0)
	eqcounts = []
	for i in equal:
		classification[(i,i)] = IntronType.equal
		eqcounts.append(counts[(i[0]-1,i[1])])
	eq = len(equal)
	extra = short.difference(chosen)
	excounts = []
	for i in extra:
		excounts.append(counts[(i[0]-1,i[1])])
	
	ex_only = len(extra.intersection(ref))		
	#print("Counts in equal introns: ", eqcounts)
	#print("Counts in extra introns: ", excounts)
	#print("Counts in paired introns: ", scounts)
	return classification, eq, excounts, scounts, ex_only


start = time.time()
args = parse_args()
#print(args.short)
short_int = get_introns(args.short)
#print(short_int[1])
if(args.iso[-2:] != "db"):
	gtf = args.iso
	args.output_exists = os.path.exists(args.out)
	if not args.output_exists:
		os.makedirs(args.out)
	db_file = os.path.join(args.out, "mouse_isoquant.db")
	gtf2db(gtf, db_file, True)
else:
	db_file = args.iso
if(args.ref[-2:] != "db"):
	gtf = args.ref
	args.output_exists = os.path.exists(args.out)
	if not args.output_exists:
		os.makedirs(args.out)
	ref = os.path.join(args.out, "mouse_reference.db")
	gtf2db(gtf, ref, True)
else:
	ref = args.ref
db = gffutils.FeatureDB(db_file)
ref_db = gffutils.FeatureDB(ref)
reference = introns_from_db(ref_db)
gene_list = list(db.features_of_type('gene', order_by=('seqid', 'start')))
cl = {}
onlyshort = 0
current_region = (0,0)
chromosome = 0
genes = []
equal_c = 0
extra_c = []
paired_c = []
for gene in gene_list:
	if overlaps((gene.start,gene.end), current_region) and gene.seqid == chromosome: # Chromosom mit abfragen?
		genes.append(gene)
		current_region = (current_region[0], max(current_region[1], gene.end))
	else:
		if genes : 
			iso_int, short_int, counts = introns_from_region(db, genes, current_region, args.short)
			ref_int = ref_from_region(ref_db, (chromosome, current_region[0], current_region[1]))
			#print(len(iso_int.difference(ref_int)))
			c, eq, ex, p, s = compare_alg(short_int, iso_int, counts, ref_int, chromosome)
			cl.update(c)
			#equal_c.extend(eq)
			equal_c += eq
			extra_c.extend(ex)
			paired_c.extend(p)
			onlyshort += s
		genes = [gene]
		current_region = (gene.start, gene.end)
		chromosome = gene.seqid
# final dataset
print("Number of correct introns only found in Illumina:", onlyshort)
print("Number of Introns equal in IsoQuant and Illumina:", equal_c)
end = time.time()
print(end - start)

