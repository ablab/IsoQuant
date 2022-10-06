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

BAMBU_GTF=""

TALON_GTF=""

FLAIR_GTF=""

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
