#!/bin/bash
# ─────────────────────────────────────────────────────────────
# 01_align.sh — Phase 1: Alignment Pipeline
# C. elegans paired-end ChIP-seq
#
# Usage: bash 01_align.sh <PE1.fastq.gz> <PE2.fastq.gz> <sample_name>
# Example: bash 01_align.sh rep1_PE1.fastq.gz rep1_PE2.fastq.gz rep1
# ─────────────────────────────────────────────────────────────

PE1=$1
PE2=$2
SAMPLE=$3

# ── Edit this path to point to your bowtie2 index prefix ──────
# Libuda Lab: /projects/libudalab/jbrilon/N2_genome_index
GENOME_INDEX="/projects/libudalab/jbrilon/N2_genome_index"

ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
QUALITY=20
MIN_LENGTH=20
MAPQ=30
THREADS=8

# ── cutadapt ──────────────────────────────────────────────────
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
picard MarkDuplicates \
  I=${SAMPLE}_sorted.bam \
  O=${SAMPLE}_deduped.bam \
  M=${SAMPLE}_dup_metrics.txt \
  REMOVE_DUPLICATES=true

samtools index ${SAMPLE}_deduped.bam

# ── sambamba ──────────────────────────────────────────────────
sambamba view \
  -h \
  -f bam \
  -F "mapping_quality >= $MAPQ and not unmapped and not (ref_name == 'chrM') and proper_pair" \
  -o ${SAMPLE}_filtered.bam \
  ${SAMPLE}_deduped.bam

samtools index ${SAMPLE}_filtered.bam
