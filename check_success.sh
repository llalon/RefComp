#!/bin/bash

# This script checks for completion of each step.

source env.sh

# Check if jobs are done
for file in $ROOT_DIR/*.out; do
    [ -e "$file" ] || continue
    
    if [[ 'grep "All done now" $file' ]]; then
        # Some Actions
        echo 'job '$file' done'
    else
        echo 'job '$file' not done'
    fi
    
done

# check if all 40 alignments went through
if [[ $(ls -1q $OUT_DIR/*.sam | wc -l) = "40" ]]; then
    echo "All 40 alignments completed"
fi


# check if all 40 stats went through (40x2 ea)
if [[ $(ls -1q $STAT_DIR/*.tsv | wc -l) = "80" ]]; then
    echo "Stats calculated on all 40 alignments"
fi