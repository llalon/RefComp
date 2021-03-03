#!/bin/sh
#
# This script will run all the required sbatch scripts

# Source our variables
source env.sh

# Make sure folders exist
mkdir -p $OUT_DIR

# Only index if not already complete
index=true

if [ "$index" = false ]; then
    echo "Index not complete. Adding job to queue. Please run this script with the option index=true once complete"
    
    ## Step. 1 - index
    for ref in $REF_BUR $REF_COD; do
        [ -e "$ref" ] || continue
        
        # Index for both bwa and bowtie2
        echo "ADDING TO QUEUE: bwa and bowtie2 indexing of $ref"
        sbatch $SCRQ_DIR/run_bwabowtie2index_queuesub.sh $ref
        
    done
    
else
    
    ## Step. 2 - Run alignemnts. 10x per reference.
    for fastq in $GEN_MATCH/*.fq.gz; do
        [ -e "$fastq" ] || continue
        
        for ref in $REF_BUR $REF_COD; do
            [ -e "$ref" ] || continue
            
            # bwa
            echo "ADDING TO QUEUE: bwa alignment of $fastq to $ref"
            sbatch $SCRQ_DIR/run_bwa_queuesub.sh $ref $fastq $OUT_DIR
            
            # bowtie2
            echo "ADDING TO QUEUE: bowtie2 alignment of $fastq to $ref"
            sbatch $SCRQ_DIR/run_bowtie2_queuesub.sh $ref $fastq $OUT_DIR
            
        done
    done
    
fi

