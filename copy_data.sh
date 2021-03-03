#!/bin/bash
#
# copy data file to this dir to keep everything in same folder
#
# 10x individual burbot raw sequence data
# 1x burbot reference genome (low quality, fragmented)
# 1x cod reference genome (high quality)

source env.sh

mkdir -p $DATA_DIR

cp -rv /scratch/emandevi/genomic_methods_w2021/Project2/burbot_raw_data $DATA_DIR
cp -rv /scratch/emandevi/genomic_methods_w2021/Project2/burbot_reference_genome $DATA_DIR
cp -rv /scratch/emandevi/genomic_methods_w2021/Project2/cod_reference_genome $DATA_DIR
