#!/usr/bin/env bash
# dataset settings
DATASET_NAME="ONT_R10"
DATA_DIR="data/"
REF_GENOME=$DATA_DIR"GRCm39.primary_assembly.genome.fa"
REDUCED_ANNOTATION_PREFIX=$DATA_DIR"mouse.gencode.M26.spatial.20percent"
REDUCED_ANNOTATION=$REDUCED_ANNOTATION_PREFIX".reduced.gtf"
REF_ANNOTATION=$REDUCED_ANNOTATION
INPUT_READS=$DATA_DIR"Mouse.ONT.R10.4.simulated.fasta"
INPUT_BAM=$DATA_DIR"Mouse.ONT.R10.4.simulated.sorted.bam"
ISOQUANT_DATA_TYPE="nanopore"
THREADS="20"

# paths to tools, end with slash
ISOQUANT_PATH=""
STRINGTIE_PATH=""
FLAIR_PATH=""
BAMBU_PATH=""
TALON_PATH=""

# BAM file can also be obtained with the following command
# minimap2 must be in $PATH
if [ -z "$INPUT_BAM" ]
then
  $ISOQUANT_PATH/isoquant.py --reference $REF_GENOME --complete_genedb --genedb $REF_ANNOTATION --fastq $INPUT_READS \
  -l $DATASET_NAME -d $ISOQUANT_DATA_TYPE --run_aligner_only -t $THREADS -o "mapped_reads_"$DATASET_NAME
  INPUT_BAM=`realpath "mapped_reads_"$DATASET_NAME/$DATASET_NAME/aux/*.bam`
fi

# reference based transcript discovery

$ISOQUANT_PATH"isoquant.py" --reference $REF_GENOME --complete_genedb --genedb $REF_ANNOTATION --bam $INPUT_BAM \
-d $ISOQUANT_DATA_TYPE -t $THREADS -l $DATASET_NAME -o "IsoQuant_out_"$DATASET_NAME
ISOQUANT_GTF="IsoQuant_out_"$DATASET_NAME/$DATASET_NAME/$DATASET_NAME".transcript_models.gtf"

STRINGTIE_GTF="StringTie_out_"$DATASET_NAME".gtf"
$STRINGTIE_PATH"stringtie" -G $REF_ANNOTATION -L $INPUT_BAM -p $THREADS -o $STRINGTIE_OUT

BAMBU_GTF="Bambu_out_"$DATASET_NAME".gtf"
awk ' $3 >= 1 ' counts_transcript.txt | sort -k3,3n > expressed_annotations.gtf.counts
cut -f1 expressed_annotations.gtf.counts > expressed_transcripts.txt
grep -Ff expressed_transcripts.txt extended_annotations.gtf > $BAMBU_GTF

FLAIR_PREFIX="Flair_out_"$DATASET_NAME
FLAIR_GTF=$FLAIR_PREFIX".gtf"
$FLAIR_PATH"bin/bam2Bed12.py" -i $INPUT_BAM > $FLAIR_PREFIX".bed"
$FLAIR_PATH"flair.py" correct -q $FLAIR_PREFIX".bed" -g $REF_GENOME -f $REF_ANNOTATION -o $FLAIR_PREFIX -t $THREADS
$FLAIR_PATH"flair.py" collapse -g $REF_GENOME -f $REF_ANNOTATION -r $INPUT_READS -q $FLAIR_PREFIX"_all_corrected.bed" -o $FLAIR_PREFIX -t $THREADS

#TALON requires config file: dataset name, sample description, platform, sam file (comma-delimited)
TALON_PREFIX="Talon_out_"$DATASET_NAME
TALON_GTF=$TALON_PREFIX".gtf"
$TALON_PATH"talon_label_reads" --f $INPUT_BAM --t=$THREADS --o=$TALON_PREFIX --g $REF_GENOME
samtools calmd -@ $THREADS --reference $REF_GENOME $TALON_PREFIX"_labeled.sam"  > $TALON_PREFIX"_labeled.md.sam"
$TALON_PATH"talon_initialize_database" --f $REF_ANNOTATION --g $DATASET_NAME --a $DATASET_NAME --o $DATASET_NAME
echo $DATASET_NAME","$DATASET_NAME",ONT,"$TALON_PREFIX"_labeled.md.sam" > $TALON_PREFIX".csv"
$TALON_PATH"talon" --build $DATASET_NAME --db $DATASET_NAME".db" --o $TALON_PREFIX"_raw" --f $TALON_PREFIX".csv"
$TALON_PATH"talon_filter_transcripts" --db $DATASET_NAME".db" -a $DATASET_NAME --datasets $DATASET_NAME --o $TALON_PREFIX"_filter" --f $TALON_CSV
$TALON_PATH"talon_create_GTF" --build $DATASET_NAME --db $DATASET_NAME".db" -a $DATASET_NAME --o $TALON_PREFIX  --whitelist=$TALON_PREFIX"_filter"

# reduced annotation analysis, path to gffcompare must be in $PATH
$ISOQUANT_PATH/misc/reduced_db_gffcompare.py --genedb $REDUCED_ANNOTATION_PREFIX --gtf $ISOQUANT_GTF --tool isoquant -o $DATASET_NAME"_isoquant_reduced_db"
$ISOQUANT_PATH/misc/reduced_db_gffcompare.py --genedb $REDUCED_ANNOTATION_PREFIX --gtf $STRINGTIE_GTF --tool stringtie -o $DATASET_NAME"_stringtie_reduced_db"
# for Bambu counts must be located in $BAMBU_GTF.counts
$ISOQUANT_PATH/misc/reduced_db_gffcompare.py --genedb $REDUCED_ANNOTATION_PREFIX --gtf $BAMBU_GTF --tool bambu -o $DATASET_NAME"_bambu_reduced_db"
$ISOQUANT_PATH/misc/reduced_db_gffcompare.py --genedb $REDUCED_ANNOTATION_PREFIX --gtf $TALON_GTF --tool talon -o $DATASET_NAME"_talon_reduced_db"
$ISOQUANT_PATH/misc/reduced_db_gffcompare.py --genedb $REDUCED_ANNOTATION_PREFIX --gtf $FLAIR_GTF --tool falir -o $DATASET_NAME"_flair_reduced_db"


# de novo annotation analysis, path to gffcompare must be in $PATH
echo $TALON_GTF >> gtfs.list
echo $FLAIR_GTF >> gtfs.list
echo $BAMBU_GTF >> gtfs.list
echo $STRINGTIE_GTF >> gtfs.list
echo $ISOQUANT_PATH >> gtfs.list
$ISOQUANT_PATH/misc/denovo_model_stats.py --gtf_list gtfs.list -o $DATASET_NAME"_denovo_stats"

# reference-free analysis, path to gffcompare must be in $PATH
$ISOQUANT_PATH"isoquant.py" --reference $REF_GENOME --bam $INPUT_BAM \
-d $ISOQUANT_DATA_TYPE -t $THREADS -l $DATASET_NAME -o "IsoQuant_denovo_out_"$DATASET_NAME
ISOQUANT_GTF="IsoQuant_denovo_out_"$DATASET_NAME/$DATASET_NAME/$DATASET_NAME".transcript_models.gtf"

STRINGTIE_GTF="StringTie_denovo_out_"$DATASET_NAME".gtf"
$STRINGTIE_PATH"stringtie" -L $INPUT_BAM -p $THREADS -o $STRINGTIE_OUT

gffcompare -r $REDUCED_ANNOTATION_PREFIX".expressed.gtf" -o  $DATASET_NAME"_isoquant_gffcompare" ISOQUANT_GTF
gffcompare -r $REDUCED_ANNOTATION_PREFIX".expressed.gtf" -o  $DATASET_NAME"_stringtie_gffcompare" $STRINGTIE_GTF
