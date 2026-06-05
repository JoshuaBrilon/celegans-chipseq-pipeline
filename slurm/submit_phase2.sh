#!/bin/bash
#SBATCH --job-name=chipseq_phase2
#SBATCH --partition=computelong
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/phase2_%j.out
#SBATCH --error=logs/phase2_%j.err

module load samtools/1.19

# Update this path to your conda installation
export PATH=/home/${USER}/miniconda3/envs/chipseq/bin:$PATH

# Update this path to your data directory
cd /projects/libudalab/${USER}/

mkdir -p logs signal peaks idr
bash scripts/02_phase2.sh
