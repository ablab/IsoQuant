#!/bin/bash
############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# Generate the partial runs the RESUME* configs resume from (see .claude/TESTING_SYSTEM.md).
# Every run is killed, with all of its worker processes, at a chosen log line (run_until.py);
# the script fails if a run ends before its stop point, so a finished run is never stored
# as a partial one. Regenerate after any change to the checkpoint layout
# (.claude/RESUME_CHECKPOINTS.md): checkpoints of other versions are refused with exit code 26.
#
# Usage: generate_resume_test.sh <isoquant_dir> [sirv|sc|all] [resume_data_dir]
# The SC data comes from prepare_sc_resume_data.sh.

set -euo pipefail

ISOQUANT_DIR=$(realpath "${1:?isoquant directory}")
WHICH=${2:-all}
RESUME_DIR=${3:-/abga/work/andreyp/ci_isoquant/data/resume_checkpoints}
DATA=/abga/work/andreyp/ci_isoquant/data
RUN_UNTIL="python3 $ISOQUANT_DIR/isoquant_tests/github/run_until.py"
ISOQUANT="python3 $ISOQUANT_DIR/isoquant.py"
mkdir -p "$RESUME_DIR"

# partial <name> <run_until stop options...> -- <isoquant options...>
partial() {
    local name=$1; shift
    local stop_opts=()
    while [ "$1" != "--" ]; do stop_opts+=("$1"); shift; done
    shift
    rm -rf "${RESUME_DIR:?}/$name"
    echo "=== $name"
    $RUN_UNTIL "${stop_opts[@]}" --log "$RESUME_DIR/$name.run_until.log" -- \
        $ISOQUANT -o "$RESUME_DIR/$name" "$@" > /dev/null
}

# complete <name> <isoquant options...>: a finished run, resuming it must change nothing
complete() {
    local name=$1; shift
    rm -rf "${RESUME_DIR:?}/$name"
    echo "=== $name"
    $ISOQUANT -o "$RESUME_DIR/$name" "$@" > "$RESUME_DIR/$name.run.log" 2>&1
}

generate_sirv() {
    local opts=(-r $DATA/ref/SIRV4/SIRV_isoforms_multi-fasta_200709a.fasta -d ont -p SRIV.Set4.ONT_R10
                --genedb $DATA/ref/SIRV4/SIRV_isoforms_multi-fasta-annotation.reduced.gtf
                -t 6 --complete_genedb --force)
    local bam=/abga/work/andreyp/data/sirvs/Lexogen.SIRVs.Set4.ONT_cDNA.R10.4.sorted.bam
    # in read collection, one chromosome done
    partial SRIV.Set4.ONT_R10.resume1 --after "Collecting read alignments" --stop "Finished processing chromosome" \
        -- "${opts[@]}" --bam $bam
    # in model construction
    partial SRIV.Set4.ONT_R10.resume2 --stop "Gene and transcript records have unequal strands: SIRV5: -, SIRV504: +" \
        -- "${opts[@]}" --bam $bam
    # read collection finished (marker written, before construction)
    partial SRIV.Set4.ONT_R10.resume3 --stop "To keep these intermediate files for" \
        -- "${opts[@]}" --bam $bam
    # in model construction, with extra outputs
    partial SRIV.Set4.ONT_R10.resume4 --after "Processing assigned reads" --stop "Finished processing chromosome" \
        -- "${opts[@]}" --bam $bam --count_exons --sqanti_output --check_canonical
    # in read collection, two BAM files
    partial SRIV.Set4.ONT_R10.resume_group1 --after "Collecting read alignments" --stop "Finished processing chromosome" \
        -- "${opts[@]}" --bam $DATA/groups/Lexogen.SIRVs.Set4.ONT_cDNA.R10.4.grouped.group0.bam \
        $DATA/groups/Lexogen.SIRVs.Set4.ONT_cDNA.R10.4.grouped.group1.bam
    # in read collection, read groups from a table
    partial SRIV.Set4.ONT_R10.resume_group2 --after "Collecting read alignments" --stop "Finished processing chromosome" \
        -- "${opts[@]}" --bam $DATA/groups/Lexogen.SIRVs.Set4.ONT_cDNA.R10.4.sorted.bam \
        --read_group file:$DATA/groups/Lexogen.SIRVs.Set4.ONT_cDNA.R10.4.grouped.groups.tsv
}

generate_sc() {
    local sc=$DATA/sc_resume
    local opts=(-r $sc/Mouse.10x.4chr.genome.fa --genedb $sc/Mouse.10x.4chr.gtf -d ont
                -m tenX_v3 --split_molecules false --complete_genedb
                --barcode_whitelist $DATA/barcodes/10xMultiome_5K.tsv
                --large_output allinfo read2transcripts tagged_bam deduplicated_bam -t 4 --force)
    local one=(-p Mouse.10x.4chr --bam $sc/Mouse.10x.4chr.bam)
    # in barcode calling, nothing marked yet
    partial Mouse.10x.4chr.resume_sc1 --stop "Minimal alignment score set to" -- "${opts[@]}" "${one[@]}"
    # barcodes called and split, in the tagged BAM (2 of 4 fragments written)
    partial Mouse.10x.4chr.resume_sc2 --after "Writing tagged BAM" --stop "alignments for chromosome" --occurrence 2 \
        -- "${opts[@]}" "${one[@]}"
    # in read collection, 2 of 4 chromosomes done
    partial Mouse.10x.4chr.resume_sc3 --after "Collecting read alignments" --stop "Finished processing chromosome" \
        --occurrence 2 -- "${opts[@]}" "${one[@]}"
    # all chromosomes collected, multimappers not resolved, collect not marked
    partial Mouse.10x.4chr.resume_sc4 --stop "Resolving multimappers" -- "${opts[@]}" "${one[@]}"
    # in UMI filtering, 2 of 4 chromosomes done
    partial Mouse.10x.4chr.resume_sc5 --stop "PCR duplicates filtered for chromosome" --occurrence 2 \
        -- "${opts[@]}" "${one[@]}"
    # UMI filtering done, in the deduplicated BAM
    partial Mouse.10x.4chr.resume_sc6 --after "Writing deduplicated BAM" --stop "Deduplicated" --occurrence 2 \
        -- "${opts[@]}" "${one[@]}"
    # in assignment processing / construction, 2 of 4 chromosomes done
    partial Mouse.10x.4chr.resume_sc7 --after "Processing assigned reads" --stop "Finished processing chromosome" \
        --occurrence 2 -- "${opts[@]}" "${one[@]}"
    # in the merge: gene counts converted, transcript counts not
    partial Mouse.10x.4chr.resume_sc8 --stop "Matrix was saved to 3 files" -- "${opts[@]}" "${one[@]}"
    # all outputs and the summary written, experiment not marked complete, nothing cleaned up
    partial Mouse.10x.4chr.resume_sc9 --stop "Run summary is stored in" -- "${opts[@]}" "${one[@]}"
    # a finished run
    complete Mouse.10x.4chr.resume_sc10 "${opts[@]}" "${one[@]}"
    # two experiments: S1 complete, S2 in read collection
    partial Mouse.10x.4chr.resume_sc_2s --after "Processing experiment S2" --stop "Finished processing chromosome" \
        -- "${opts[@]}" --yaml $sc/Mouse.10x.4chr.2samples.yaml
}

case "$WHICH" in
    sirv) generate_sirv ;;
    sc) generate_sc ;;
    all) generate_sirv; generate_sc ;;
    *) echo "unknown set: $WHICH (sirv|sc|all)"; exit 2 ;;
esac
echo "Partial runs are in $RESUME_DIR"
