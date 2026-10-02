#!/bin/bash
############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# One-time preparation of the single-cell resume test data (RESUME_SC* configs):
# ~100k reads of the 10x mouse dataset on four small chromosomes, the matching genome
# and annotation, and a 2-experiment split of the same reads for the YAML variant.
#
# Usage: prepare_sc_resume_data.sh [output_dir]

set -euo pipefail

OUT=${1:-/abga/work/andreyp/ci_isoquant/data/sc_resume}
SRC_BAM=/abga/work/andreyp/ci_isoquant/data/barcodes/Mouse.10x.5k.ONT_cDNA.R10.4.no_trunc.bam
SRC_GENOME=/abga/work/andreyp/data/reference/mouse/GRCm39.primary_assembly.genome.fa
SRC_GTF=/abga/work/andreyp/data/reference/mouse/gencode.vM36.basic.annotation.gtf
CHROMOSOMES="chr14 chr16 chr18 chr19"
N_READS=100000
SEED=42

mkdir -p "$OUT"
cd "$OUT"
PREFIX=Mouse.10x.4chr

# genome and annotation restricted to the four chromosomes
samtools faidx "$SRC_GENOME" $CHROMOSOMES > $PREFIX.genome.fa
samtools faidx $PREFIX.genome.fa
awk -v chrs="$CHROMOSOMES" 'BEGIN { n = split(chrs, a, " "); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
     /^#/ || ($1 in keep)' "$SRC_GTF" > $PREFIX.gtf

# a deterministic random subset of read names with an alignment on these chromosomes
samtools view -F 0x900 "$SRC_BAM" $CHROMOSOMES | cut -f1 | sort -u > all_reads.txt
python3 - "$N_READS" "$SEED" <<'EOF'
import random, sys
n, seed = int(sys.argv[1]), int(sys.argv[2])
names = [l.rstrip("\n") for l in open("all_reads.txt")]
random.Random(seed).shuffle(names)
picked = sorted(names[:n])
with open("reads.txt", "w") as f:
    f.write("\n".join(picked) + "\n")
# the 2-experiment variant: same reads, split in halves
half = len(picked) // 2
with open("reads.S1.txt", "w") as f:
    f.write("\n".join(picked[:half]) + "\n")
with open("reads.S2.txt", "w") as f:
    f.write("\n".join(picked[half:]) + "\n")
EOF

# every alignment of a picked read on these chromosomes, header reduced to them
samtools view -H "$SRC_BAM" | awk -v chrs="$CHROMOSOMES" '
    BEGIN { n = split(chrs, a, " "); for (i = 1; i <= n; i++) keep["SN:" a[i]] = 1 }
    $1 != "@SQ" || ($2 in keep)' > header.sam
subset_bam() {
    local names=$1 out=$2
    samtools view -N "$names" "$SRC_BAM" $CHROMOSOMES | cat header.sam - | samtools sort -o "$out" -
    samtools index "$out"
}
subset_bam reads.txt $PREFIX.bam
subset_bam reads.S1.txt $PREFIX.S1.bam
subset_bam reads.S2.txt $PREFIX.S2.bam

cat > $PREFIX.2samples.yaml <<EOF
[
  data format: "bam",
  {
    name: "S1",
    long read files: ["$OUT/$PREFIX.S1.bam"]
  },
  {
    name: "S2",
    long read files: ["$OUT/$PREFIX.S2.bam"]
  }
]
EOF

rm -f all_reads.txt header.sam
echo "Reads: $(wc -l < reads.txt), experiment S1: $(wc -l < reads.S1.txt), S2: $(wc -l < reads.S2.txt)"
samtools idxstats $PREFIX.bam
