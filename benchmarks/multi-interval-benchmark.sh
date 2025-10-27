#!/bin/bash

#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

module load java
module load gatk
module load picard
module load samtools
module load sambamba

# Compare the runtime of different downsampling approach from various tools
# Assumes GATK, samtools, and sambamba are installed and available
# Assumes Python is available and samsampleX hass all dependencies installed in a virtual environment

# This script will run the downsampling approach for each tool and print the runtime to a file
# The chromosomal coordinates can be adjusted to the region of interest

# Default coordinates are a padded range for the Human Leukocyte Antigen (HLA) region
# By default, the region is divided into intervals of 1000 bases

CHR=chr6
START=25000000
END=35000000
INTERVAL_LENGTH=1000
INPUT_BAM=HG002.GRCh38.300x.bam
SEED=42

INPUTS=benchmarks/gatk-inputs
OUTPUTS=benchmarks/outputs
mkdir -p $INPUTS
mkdir -p $OUTPUTS

# GATK downsampling cannot limit the operation to a specific interval, operating on the entire BAM file
# The input BAM file is thus pre-processed using samtools first, which is not included in benchmarking
echo "Pre-processing BAM files for GATK..."
for i in $(seq ${START} ${INTERVAL_LENGTH} $((END-INTERVAL_LENGTH)));
do
    INTERVAL_BAM=$INPUTS/${CHR}-${i}-$((i+INTERVAL_LENGTH)).bam
    samtools view $INPUT_BAM $CHR:$i-$((i+INTERVAL_LENGTH)) -b -o $INTERVAL_BAM && \
    samtools index $INTERVAL_BAM
done

echo "Benchmarking GATK HighAccuracy..."
GATK_HIGH_ACCURACY=$OUTPUTS/multi-interval-gatk-high-accuracy.log
env time -v -o $GATK_HIGH_ACCURACY ./benchmarks/gatk-highaccuracy.sh $CHR $START $END $INTERVAL_LENGTH $INPUT_BAM $SEED $INPUTS $OUTPUTS

echo "Benchmarking GATK ConstantMemory..."
GATK_CONSTANT_MEMORY=$OUTPUTS/multi-interval-gatk-constant-memory.log
env time -v -o $GATK_CONSTANT_MEMORY ./benchmarks/gatk-constant-memory.sh $CHR $START $END $INTERVAL_LENGTH $INPUT_BAM $SEED $INPUTS $OUTPUTS

echo "Benchmarking GATK Chained..."
GATK_CHAINED=$OUTPUTS/multi-interval-gatk-chained.log
env time -v -o $GATK_CHAINED ./benchmarks/gatk-chained.sh $CHR $START $END $INTERVAL_LENGTH $INPUT_BAM $SEED $INPUTS $OUTPUTS

# for i in $(seq ${START} ${INTERVAL_LENGTH} $((END-INTERVAL_LENGTH)));
# do
#     gatk AddOrReplaceReadGroups -I $INPUTS/$CHR-$i-$((i+INTERVAL_LENGTH)).bam -O $INPUTS/$CHR-$i-$((i+INTERVAL_LENGTH))-with-read-groups.bam -LB 1 -PL ILLUMINA -PU unit1 -SM HG002
# done

# echo "Benchmarking GATK DownsampleByDuplicateSet..."
# GATK_DOWN_BY_DUP_SET=$OUTPUTS/multi-interval-gatk-downsample-by-dup-set.log
# env time -o $GATK_DOWN_BY_DUP_SET -v bash -c "
#         for i in \$(seq \${START} \${INTERVAL_LENGTH} \$((END-INTERVAL_LENGTH))); do
#             gatk DownsampleByDuplicateSet --fraction-to-keep 0.1 -I $INPUTS/\$CHR-\$i-\$((i+INTERVAL_LENGTH))-with-read-groups.bam -O $OUTPUTS/gatk-downsample-by-dup-set.bam
#         done
#     "

echo "Benchmarking samtools..."
SAMTOOLS=$OUTPUTS/multi-interval-samtools.log
env time -v -o $SAMTOOLS ./benchmarks/samtools.sh $CHR $START $END $INTERVAL_LENGTH $INPUT_BAM $SEED $INPUTS $OUTPUTS

SAMBAMBA=$OUTPUTS/multi-interval-sambamba.log
echo "Benchmarking sambamba..."
env time -v -o $SAMBAMBA ./benchmarks/sambamba.sh $CHR $START $END $INTERVAL_LENGTH $INPUT_BAM $SEED $INPUTS $OUTPUTS

# samsampleX mapping is also done outside of the benchmarking process
# Assumes samsampleX is installed and available
samsampleX map --in-bam $INPUT_BAM --contig $CHR --start $START --end $END --interval-length $INTERVAL_LENGTH --bed $OUTPUTS/multi-interval-benchmark.bed && \

echo "Benchmarking samsampleX..."
SUBSAMPLE_READS=$OUTPUTS/multi-interval-samsampleX.log
env time -o $SUBSAMPLE_READS -v samsampleX sample --seed $SEED --in-bam $INPUT_BAM --bed $OUTPUTS/multi-interval-benchmark.bed --out-bam $OUTPUTS/samsampleX.bam

get_stat(){
    grep $1 $2 | rev | cut -d' ' -f1 | rev
}

OUTPUT=benchmarks/multi-interval-benchmark.tsv
make_table() {

    echo "Tool/Mode	User time	System time	Wall clock  Memory used" > $OUTPUT
    echo "gatk-HighAccuracy	$(get_stat User $GATK_HIGH_ACCURACY)	$(get_stat System $GATK_HIGH_ACCURACY)	$(get_stat wall $GATK_HIGH_ACCURACY)    $(get_stat Maximum $GATK_HIGH_ACCURACY)KB" >> $OUTPUT
    echo "gatk-ConstantMemory	$(get_stat User $GATK_CONSTANT_MEMORY)	$(get_stat System $GATK_CONSTANT_MEMORY)	$(get_stat wall $GATK_CONSTANT_MEMORY)  $(get_stat Maximum $GATK_CONSTANT_MEMORY)KB" >> $OUTPUT
    echo "gatk-Chained	$(get_stat User $GATK_CHAINED)	$(get_stat System $GATK_CHAINED)	$(get_stat wall $GATK_CHAINED)  $(get_stat Maximum $GATK_CHAINED)KB" >> $OUTPUT
    # echo "gatk-DownsampleByDuplicateSet	$(get_stat User $GATK_DOWN_BY_DUP_SET)	$(get_stat System $GATK_DOWN_BY_DUP_SET)	$(get_stat wall $GATK_DOWN_BY_DUP_SET)  $(get_stat Maximum $GATK_DOWN_BY_DUP_SET)KB" >> $OUTPUT
    echo "samtools	$(get_stat User $SAMTOOLS)	$(get_stat System $SAMTOOLS)	$(get_stat wall $SAMTOOLS)   $(get_stat Maximum $SAMTOOLS)KB" >> $OUTPUT
    echo "sambamba	$(get_stat User $SAMBAMBA)	$(get_stat System $SAMBAMBA)	$(get_stat wall $SAMBAMBA)  $(get_stat Maximum $SAMBAMBA)KB" >> $OUTPUT
    echo "samsampleX	$(get_stat User $SUBSAMPLE_READS)	$(get_stat System $SUBSAMPLE_READS)	$(get_stat wall $SUBSAMPLE_READS) $(get_stat Maximum $SUBSAMPLE_READS)KB" >> $OUTPUT
}

make_table
echo "Benchmarking completed"