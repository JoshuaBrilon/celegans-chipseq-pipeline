#!/bin/bash
#SBATCH --job-name=chipseq_ctrl
#SBATCH --partition=computelong
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/ctrl_%j.out
#SBATCH --error=logs/ctrl_%j.err

# ─────────────────────────────────────────────────────────────
# SLURM wrapper for control BAM creation
#
# Usage: sbatch slurm/submit_ctrl.sh <ctrl_PE1.fastq.gz> <ctrl_PE2.fastq.gz>
#
# What is this?
# Your control (input) sample needs to be processed through the
# same alignment pipeline as your ChIP samples before it can be
# used in peak calling. This script handles that automatically.
#
# After this job finishes, you will have ctrl_filtered.bam which
# is required by 02_phase2.sh.
# ─────────────────────────────────────────────────────────────

PE1=$1
PE2=$2

module load samtools/1.19
module load bowtie2/2.4.4
module load cutadapt/4.9-20240729

# Update this path to your conda installation
export PATH=/home/${USER}/miniconda3/envs/chipseq/bin:$PATH
conda activate chipseq

# Update this path to your data directory
cd /projects/libudalab/${USER}/

mkdir -p logs
bash scripts/00_create_ctrl_bam.sh $PE1 $PE2
