#!/bin/sh

# Adapted from Class Script
# USAGE: run_bwa.sh reference fasta output_dir

#SBATCH --account=def-emandevi
#SBATCH --time=0-00:15:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)

## Load modules
module load bowtie2
module load samtools/1.9

# Parse args
# Reference genome (abs path)
REF=$1
# Fastq to align (abs path)
FASTQ=$2
# Output directory to put files in (abs path)
OUT_DIR=$3

ref_basename=`echo $(basename $REF) | sed 's/\.f[a-z]*//g'`
basename=`echo $(basename $FASTQ) | sed 's/\.f[a-z]*//g'`
ref_path=$(dirname $REF)
output_file=$OUT_DIR/$basename.bowtie2.$ref_basename

echo "Aligning $REF using bowtie2..."

# Align
echo "Starting alignment of $FASTQ to $REF"
bowtie2 -x $ref_path'/'$ref_basename -U $FASTQ -S $output_file.sam --very-sensitive-local

# Convert
echo "Converting sam to bam for $basename"
samtools view -b -S -o $output_file.bam $output_file.sam

# Sort
echo "Sorting and indexing bam files for $basename"
samtools sort $output_file.bam -o $output_file.sorted.bam
samtools index $output_file.sorted.bam

# Fix/detect errors
echo "Cleaning up the mess... just a minute!"
if ![[ -s $output_file.bam ]]; then
    echo "$basename.sam is empty! Something's fishy..."
fi

if ![[ -s $output_file.bam ]]; then
    echo "$basename.bam is empty! Something's fishy..."
fi

echo "All done now"