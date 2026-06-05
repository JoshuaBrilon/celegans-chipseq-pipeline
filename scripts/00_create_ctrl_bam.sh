#!/bin/bash
# ─────────────────────────────────────────────────────────────
# 00_create_ctrl_bam.sh — Create control (input) BAM for ChIP-seq
#
# What is a control BAM?
# ─────────────────────────────────────────────────────────────
# In a ChIP-seq experiment, the "control" or "input" is a sample
# prepared identically to your ChIP sample, but WITHOUT the
# immunoprecipitation step (no antibody pulldown). This captures
# background chromatin accessibility across the genome.
#
# MACS2 uses the control to distinguish real transcription factor
# binding sites from regions that are just generally open/accessible.
# Without a control, you cannot reliably call peaks.
#
# The control FASTQs come from the same experiment — they are the
# "input DNA" sample. Ask your PI or sequencing core which files
# are your input/control FASTQs.
#
# IMPORTANT: The control must be aligned to the SAME genome as
# your ChIP samples. This script uses the same pipeline as
# 01_align.sh to ensure that.
#
# Usage: bash 00_create_ctrl_bam.sh <PE1.fastq.gz> <PE2.fastq.gz>
# Example: bash 00_create_ctrl_bam.sh input_PE1.fastq.gz input_PE2.fastq.gz
#
# Output: ctrl_filtered.bam (used as CTRL in 02_phase2.sh)
# ─────────────────────────────────────────────────────────────

PE1=$1
PE2=$2
SAMPLE="ctrl"

# ── Edit this path to point to your bowtie2 index prefix ──────
# Libuda Lab: /projects/libudalab/jbrilon/N2_genome_index
GENOME_INDEX="/path/to/your/genome_index"

ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
QUALITY=20
MIN_LENGTH=20
MAPQ=30
THREADS=8

echo "Creating control BAM from input FASTQs..."
echo "PE1: $PE1"
echo "PE2: $PE2"
echo "Genome index: $GENOME_INDEX"
echo ""

# ── cutadapt ──────────────────────────────────────────────────
echo "Step 1/4: Trimming adapters..."
cutadapt \
  -a $ADAPTER \
  -A $ADAPTER \
  -q $QUALITY \
  -m $MIN_LENGTH \
  -o ${SAMPLE}_PE1_trimmed.fastq.gz \
  -p ${SAMPLE}_PE2_trimmed.fastq.gz \
  $PE1 \
  $PE2 \
  > ${SAMPLE}_cutadapt.log 2>&1

# ── bowtie2 + samtools ────────────────────────────────────────
echo "Step 2/4: Aligning to genome..."
bowtie2 \
  -x $GENOME_INDEX \
  -1 ${SAMPLE}_PE1_trimmed.fastq.gz \
  -2 ${SAMPLE}_PE2_trimmed.fastq.gz \
  -p $THREADS \
  --rg-id $SAMPLE \
  --rg "SM:$SAMPLE" \
  --rg "PL:ILLUMINA" \
  --rg "LB:lib1" \
  2> ${SAMPLE}_bowtie2.log \
  | samtools view -bS - \
  | samtools sort -o ${SAMPLE}_sorted.bam

samtools index ${SAMPLE}_sorted.bam

# ── picard ────────────────────────────────────────────────────
echo "Step 3/4: Removing duplicates..."
picard MarkDuplicates \
  I=${SAMPLE}_sorted.bam \
  O=${SAMPLE}_deduped.bam \
  M=${SAMPLE}_dup_metrics.txt \
  REMOVE_DUPLICATES=true

samtools index ${SAMPLE}_deduped.bam

# ── sambamba ──────────────────────────────────────────────────
echo "Step 4/4: Filtering reads..."
sambamba view \
  -h \
  -f bam \
  -F "mapping_quality >= $MAPQ and not unmapped and not (ref_name == 'chrM') and proper_pair" \
  -o ${SAMPLE}_filtered.bam \
  ${SAMPLE}_deduped.bam

samtools index ${SAMPLE}_filtered.bam

echo ""
echo "Done! Control BAM created: ctrl_filtered.bam"
echo "Use this file as CTRL in 02_phase2.sh"
echo ""
echo "Alignment summary:"
cat ${SAMPLE}_bowtie2.log
