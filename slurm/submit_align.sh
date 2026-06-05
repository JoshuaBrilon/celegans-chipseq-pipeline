#!/bin/bash
#SBATCH --job-name=chipseq_align
#SBATCH --partition=computelong
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/align_%j.out
#SBATCH --error=logs/align_%j.err

# Usage: sbatch slurm/submit_align.sh <PE1.fastq.gz> <PE2.fastq.gz> <sample_name>

PE1=$1
PE2=$2
SAMPLE=$3

module load samtools/1.19
module load bowtie2/2.4.4
module load cutadapt/4.9-20240729

# Update this path to your conda environment
export PATH=/home/${USER}/miniconda3/envs/chipseq/bin:$PATH

# Update this path to your data directory
cd /projects/libudalab/${USER}/

mkdir -p logs
bash scripts/01_align.sh $PE1 $PE2 $SAMPLE
