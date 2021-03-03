#!/bin/sh

# This script indexes a reference genome for both bwa and bowtie2. Adpated from class script.
# Usage: sbatch run_bwabowtie2index_queuesub.sh genomename.fasta

#SBATCH --account=def-emandevi
#SBATCH --time=0-03:00:00
#SBATCH --nodes=1
#SBATCH --mem=16000
#SBATCH --ntasks-per-node=1

# Load modules
module load bwa
module load bowtie2

# Creates filenames and dirs
basename=`echo $(basename $1) | sed 's/\.f[a-z]*//g'`
work_dir=$(dirname $1)

echo "Indexing $1 using bwa/bowtie2..."

# cd to the folder where the files are
cd $(dirname $work_dir)

# bwa index
echo "Running indexing for $1 with bwa"
bwa index -a bwtsw $1
echo "Done with bwa."

# bowtie2 index
echo "Indexing for bowtie2."
bowtie2-build $1 $basename
echo "Done with bowtie2."

echo "All done now"


