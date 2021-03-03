#!/bin/sh
#
# This script contains global variables for the project

# Root project dir
#ROOT_DIR='/scratch/llalon02/project_2/toastedfrog'
ROOT_DIR='/home/llalon02/projects/def-emandevi/llalon02/project2/toastedfrog'

# Data dir
DATA_DIR=$ROOT_DIR'/data'

# Directory for output files
OUT_DIR=$DATA_DIR'/output'

# Scripts dir
SCRQ_DIR=$ROOT_DIR'/sbatch_scripts'

# Genome target dir
GEN_MATCH=$DATA_DIR'/burbot_raw_data'

# Reference genomes
REF_BUR=$DATA_DIR'/burbot_reference_genome/GCA_900302385.1_ASM90030238v1_genomic.fna'
REF_COD=$DATA_DIR'/cod_reference_genome/GCF_902167405.1_gadMor3.0_genomic.fna'

# Folder to keep colates stats in
STAT_DIR=$ROOT_DIR'/stats'
