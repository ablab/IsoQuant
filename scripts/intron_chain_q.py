import sys
import pysam
from collections import defaultdict

bam_file = sys.argv[1]
quality_file = sys.argv[2]

qdict = {}
cur_id = ""
for l in open(quality_file):
    if l.startswith(">"):
        cur_id = l.strip()[1:]
    else:
        qdict[cur_id] = int(float(l.strip()))
        cur_id = ""

print("Loaded quality")
intron_len_dict = defaultdict(list)
assignment_len_dict = defaultdict(lambda: defaultdict(list))
aligned_len_dict = defaultdict(list)
c = 0
for r in pysam.AlignmentFile(bam_file, "rb").fetch():
    if r.is_supplementary or r.is_secondary or r.reference_id == -1 or not r.cigartuples:
        continue

    if r.query_name not in qdict:
        continue

    c += 1
    intron_count = sum([1 if x[0] == 3 else 0 for x in r.cigartuples])
    aligned_len = abs(r.query_alignment_start- r.query_alignment_end)

    if intron_count > 0: intron_len_dict[qdict[r.query_name]].append(intron_count)
    #if r.query_name in read_assigns: assignment_len_dict[qdict[r.query_name]][intron_count].append(read_assigns[r.query_name])
    if intron_count > 0: aligned_len_dict[qdict[r.query_name]].append(aligned_len)

    if c % 100000 == 0:
        print("Processed " + str(c))


for k in sorted(intron_len_dict.keys()):
    lst = intron_len_dict[k]
    print("%d (%d): %.4f" % (k, len(lst), sum(lst) / len(lst)))
for k in sorted(aligned_len_dict.keys()):
    lst = aligned_len_dict[k]
    print("%d (%d): %.4f" % (k, len(lst), sum(lst) / len(lst))) 
 

