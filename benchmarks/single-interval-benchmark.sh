#!/bin/bash

#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

module load java
module load gatk
module load samtools
module load sambamba

# Script to compare the runtime of different downsampling approach from various tools
# Assumes GATK, samtools, and sambamba are installed and available
# Assumes Python is available and samsampleX hass all dependencies installed in a virtual environment

CHR=chr6
START=25000000
END=35000000
INPUT_BAM=HG002.GRCh38.300x.bam
SEED=42

INPUTS=benchmarks/gatk-inputs
OUTPUTS=benchmarks/outputs
mkdir -p $INPUTS
mkdir -p $OUTPUTS

# GATK downsampling cannot limit the operation to a specific interval, operating on the entire BAM file
# The input BAM file is thus pre-processed using samtools first, not included in benchmarking
echo "Pre-processing BAM files for GATK..."
INTERVAL_BAM=$INPUTS/single-interval-subset.bam
samtools view $INPUT_BAM $CHR:$START-$END -b -o $INTERVAL_BAM && \
samtools index $INTERVAL_BAM

echo "Benchmarking GATK HighAccuracy..."
GATK_HIGH_ACCURACY=$OUTPUTS/single-interval-gatk-high-accuracy.log
env time -o $GATK_HIGH_ACCURACY -v gatk DownsampleSam -P 0.1 -R $SEED -S HighAccuracy -I $INTERVAL_BAM -O $OUTPUTS/single-interval-gatk-high-accuracy.bam

echo "Benchmarking GATK ConstantMemory..."
GATK_CONSTANT_MEMORY=$OUTPUTS/single-interval-gatk-constant-memory.log
env time -o $GATK_CONSTANT_MEMORY -v gatk DownsampleSam -P 0.1 -R $SEED -S ConstantMemory -I $INTERVAL_BAM -O $OUTPUTS/single-interval-gatk-constant-memory.bam

echo "Benchmarking GATK Chained..."
GATK_CHAINED=$OUTPUTS/single-interval-gatk-chained.log
env time -o $GATK_CHAINED -v gatk DownsampleSam -P 0.1 -R $SEED -S Chained -I $INTERVAL_BAM -O $OUTPUTS/single-interval-gatk-chained.bam

# echo "Benchmarking GATK DownsampleByDuplicateSet..."
# gatk AddOrReplaceReadGroups -I $INPUTS/single-interval-subset.bam -O $INPUTS/$CHR-$START-$END-with-read-groups.bam -LB 1 -PL ILLUMINA -PU unit1 -SM HG002
# GATK_DOWN_BY_DUP_SET=$OUTPUTS/single-interval-gatk-downsample-by-dup-set.log
# env time -o $GATK_DOWN_BY_DUP_SET -v gatk DownsampleByDuplicateSet --fraction-to-keep 0.1 -I $INPUTS/$CHR-$START-$END-with-read-groups.bam -O $OUTPUTS/single-interval-gatk-downsample-by-dup-set.bam

echo "Benchmarking samtools..."
SAMTOOLS=$OUTPUTS/single-interval-samtools.log
env time -o $SAMTOOLS -v samtools view -s $SEED.1 -b $INPUT_BAM $CHR:$START-$END -b -o $OUTPUTS/single-interval-samtools.bam

echo "Benchmarking sambamba..."
SAMBAMBA=$OUTPUTS/single-interval-sambamba.log
env time -o $SAMBAMBA -v sambamba view -s $SEED.1 $INPUT_BAM $CHR:$START-$END -o $OUTPUTS/single-interval-sambamba.bam

# samsampleX mapping is also done outside of the benchmarking process
# Assumes samsampleX is installed and available
TEMPLATE=$OUTPUTS/single-interval-benchmark.bed
samsampleX map --in-bam $INPUT_BAM --contig $CHR --start $START --end $END --interval-count 1 --bed $TEMPLATE

echo "Benchmarking samsampleX..."
SUBSAMPLE_READS=$OUTPUTS/single-interval-samsampleX.log
env time -o $SUBSAMPLE_READS -v samsampleX sample --seed $SEED --in-bam $INPUT_BAM --bed $TEMPLATE --out-bam $OUTPUTS/single-interval-samsampleX.bam

# Get last value in line from grepped file
get_stat(){
    grep $1 $2 | rev | cut -d' ' -f1 | rev
}

echo "Collecting results..."
OUTPUT=benchmarks/single-interval-benchmark.tsv
make_table() {

    echo "Tool/Mode	User time	System time	Wall clock  Memory used" > $OUTPUT
    echo "GATK HighAccuracy	$(get_stat User $GATK_HIGH_ACCURACY)	$(get_stat System $GATK_HIGH_ACCURACY)	$(get_stat wall $GATK_HIGH_ACCURACY)    $(get_stat Maximum $GATK_HIGH_ACCURACY)KB" >> $OUTPUT
    echo "GATK ConstantMemory	$(get_stat User $GATK_CONSTANT_MEMORY)	$(get_stat System $GATK_CONSTANT_MEMORY)	$(get_stat wall $GATK_CONSTANT_MEMORY)  $(get_stat Maximum $GATK_CONSTANT_MEMORY)KB" >> $OUTPUT
    echo "GATK Chained	$(get_stat User $GATK_CHAINED)	$(get_stat System $GATK_CHAINED)	$(get_stat wall $GATK_CHAINED)  $(get_stat Maximum $GATK_CHAINED)KB" >> $OUTPUT
    # echo "GATK DownsampleByDuplicateSet	$(get_stat User $GATK_DOWN_BY_DUP_SET)	$(get_stat System $GATK_DOWN_BY_DUP_SET)	$(get_stat wall $GATK_DOWN_BY_DUP_SET)  $(get_stat Maximum $GATK_DOWN_BY_DUP_SET)KB" >> $OUTPUT
    echo "samtools	$(get_stat User $SAMTOOLS)	$(get_stat System $SAMTOOLS)	$(get_stat wall $SAMTOOLS)   $(get_stat Maximum $SAMTOOLS)KB" >> $OUTPUT
    echo "sambamba	$(get_stat User $SAMBAMBA)	$(get_stat System $SAMBAMBA)	$(get_stat wall $SAMBAMBA)  $(get_stat Maximum $SAMBAMBA)KB" >> $OUTPUT
    echo "samsampleX	$(get_stat User $SUBSAMPLE_READS)	$(get_stat System $SUBSAMPLE_READS)	$(get_stat wall $SUBSAMPLE_READS) $(get_stat Maximum $SUBSAMPLE_READS)KB" >> $OUTPUT
}

make_table
echo "Benchmarking completed"