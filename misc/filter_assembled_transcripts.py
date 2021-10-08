import sys
from Bio import SeqIO

raw_transcripts_file = sys.argv[1]
filt_transcripts_file = sys.argv[2]
final_transcripts_file = sys.argv[3]

filt_transcripts = set()
with open(filt_transcripts_file) as f:
    for line in f:
        filt_transcripts.add(line.strip())

final_transcripts = []
with open(raw_transcripts_file) as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        if record.id in filt_transcripts:
            final_transcripts.append(record)

SeqIO.write(final_transcripts, final_transcripts_file, "fasta")