#!/bin/bash
#
# This script runs samtools stats against all the .bam (sorted) files in the output dir.

source env.sh

# Make sure folders exist
mkdir -p $STAT_DIR

for bam in $OUT_DIR/*.sorted.bam; do
    [ -e "$bam" ] || continue
    
    echo $(basename $bam)
    
    # get stats outputs to tsv
    samtools idxstats $bam > $STAT_DIR/$(basename $bam).reads_per_chunk.tsv
    samtools stats $bam > $STAT_DIR/$(basename $bam).stats_raw.tsv
    
done