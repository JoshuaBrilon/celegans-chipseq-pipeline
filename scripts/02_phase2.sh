#!/bin/bash
# ─────────────────────────────────────────────────────────────
# 02_phase2.sh — Phase 2: Signal generation, peak calling, IDR
# C. elegans paired-end ChIP-seq
#
# Requires:
#   rep1_filtered.bam    — replicate 1 ChIP (from 01_align.sh)
#   rep2_filtered.bam    — replicate 2 ChIP (from 01_align.sh)
#   ctrl_filtered.bam    — input/control BAM (from 00_create_ctrl_bam.sh)
#
# IMPORTANT: ctrl_filtered.bam MUST be aligned to the same genome
# as rep1 and rep2. Use 00_create_ctrl_bam.sh to create it from
# your input FASTQs. Do not use a pre-processed BAM from an
# external source unless you have verified chromosome names and
# lengths match exactly.
#
# Usage: bash 02_phase2.sh
# ─────────────────────────────────────────────────────────────

REP1="rep1_filtered.bam"
REP2="rep2_filtered.bam"
CTRL="ctrl_filtered.bam"

GENOME_SIZE="9e7"    # C. elegans effective genome size
THREADS=8

mkdir -p signal peaks idr

# ── Helper: split BAM into two random pseudoreplicates ────────
# Uses samtools -s seed.fraction to randomly subsample 50% of
# reads into each pseudoreplicate. Different seeds (1, 2) ensure
# the two pseudoreplicates are independent random draws.
split_bam() {
  local input=$1
  local pr1=$2
  local pr2=$3
  samtools view -b -s 1.5 $input | samtools sort -o $pr1
  samtools view -b -s 2.5 $input | samtools sort -o $pr2
  samtools index $pr1
  samtools index $pr2
}

# ─────────────────────────────────────────────────────────────
# STEP 1 — Signal generation (rep1)
# bamCompare computes log2(ChIP/control) at each genomic position.
# This subtracts background noise and highlights real enrichment.
# Requires ctrl_filtered.bam to be aligned to the same genome.
# ─────────────────────────────────────────────────────────────
echo "Step 1: Signal generation rep1..."

bamCompare \
  -b1 $REP1 \
  -b2 $CTRL \
  -o signal/rep1_foldchange.bigwig \
  --operation log2 \
  --normalizeUsing RPKM \
  --scaleFactorsMethod None \
  -p $THREADS

bamCompare \
  -b1 $REP1 \
  -b2 $CTRL \
  -o signal/rep1_pvalue.bigwig \
  --operation log2 \
  --normalizeUsing RPKM \
  --scaleFactorsMethod None \
  --pseudocount 1 \
  -p $THREADS

# ─────────────────────────────────────────────────────────────
# STEP 2 — Signal generation (rep2)
# ─────────────────────────────────────────────────────────────
echo "Step 2: Signal generation rep2..."

bamCompare \
  -b1 $REP2 \
  -b2 $CTRL \
  -o signal/rep2_foldchange.bigwig \
  --operation log2 \
  --normalizeUsing RPKM \
  --scaleFactorsMethod None \
  -p $THREADS

bamCompare \
  -b1 $REP2 \
  -b2 $CTRL \
  -o signal/rep2_pvalue.bigwig \
  --operation log2 \
  --normalizeUsing RPKM \
  --scaleFactorsMethod None \
  --pseudocount 1 \
  -p $THREADS

# ─────────────────────────────────────────────────────────────
# STEP 3 — Pseudoreplicated IDR (rep1)
# Split rep1 into two random halves, call peaks on each,
# then run IDR to assess self-consistency of rep1.
# ─────────────────────────────────────────────────────────────
echo "Step 3: Pseudoreplicated IDR rep1..."

split_bam $REP1 rep1_pr1.bam rep1_pr2.bam

macs2 callpeak -t rep1_pr1.bam -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n rep1_pr1 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/rep1_pr1_macs2.log
macs2 callpeak -t rep1_pr2.bam -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n rep1_pr2 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/rep1_pr2_macs2.log

sort -k8,8nr peaks/rep1_pr1_peaks.narrowPeak > peaks/rep1_pr1_sorted.narrowPeak
sort -k8,8nr peaks/rep1_pr2_peaks.narrowPeak > peaks/rep1_pr2_sorted.narrowPeak

idr --samples peaks/rep1_pr1_sorted.narrowPeak peaks/rep1_pr2_sorted.narrowPeak \
  --input-file-type narrowPeak --rank p.value \
  --output-file idr/rep1_pseudorep_idr.txt --idr-threshold 0.05 \
  --plot --log-output-file idr/rep1_pseudorep_idr.log

# ─────────────────────────────────────────────────────────────
# STEP 4 — Pseudoreplicated IDR (rep2)
# ─────────────────────────────────────────────────────────────
echo "Step 4: Pseudoreplicated IDR rep2..."

split_bam $REP2 rep2_pr1.bam rep2_pr2.bam

macs2 callpeak -t rep2_pr1.bam -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n rep2_pr1 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/rep2_pr1_macs2.log
macs2 callpeak -t rep2_pr2.bam -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n rep2_pr2 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/rep2_pr2_macs2.log

sort -k8,8nr peaks/rep2_pr1_peaks.narrowPeak > peaks/rep2_pr1_sorted.narrowPeak
sort -k8,8nr peaks/rep2_pr2_peaks.narrowPeak > peaks/rep2_pr2_sorted.narrowPeak

idr --samples peaks/rep2_pr1_sorted.narrowPeak peaks/rep2_pr2_sorted.narrowPeak \
  --input-file-type narrowPeak --rank p.value \
  --output-file idr/rep2_pseudorep_idr.txt --idr-threshold 0.05 \
  --plot --log-output-file idr/rep2_pseudorep_idr.log

# ─────────────────────────────────────────────────────────────
# STEP 5 — Pooling + signal generation
# Merge rep1 and rep2 for a combined signal track
# ─────────────────────────────────────────────────────────────
echo "Step 5: Pooling + signal generation..."

samtools merge -f pooled.bam $REP1 $REP2
samtools index pooled.bam

bamCompare \
  -b1 pooled.bam \
  -b2 $CTRL \
  -o signal/pooled_foldchange.bigwig \
  --operation log2 \
  --normalizeUsing RPKM \
  --scaleFactorsMethod None \
  -p $THREADS

bamCompare \
  -b1 pooled.bam \
  -b2 $CTRL \
  -o signal/pooled_pvalue.bigwig \
  --operation log2 \
  --normalizeUsing RPKM \
  --scaleFactorsMethod None \
  --pseudocount 1 \
  -p $THREADS

# ─────────────────────────────────────────────────────────────
# STEP 6 — Pooling + pseudoreplicated IDR
# Split pooled BAM into pseudoreplicates for pooled QC
# ─────────────────────────────────────────────────────────────
echo "Step 6: Pooled pseudoreplicated IDR..."

split_bam pooled.bam pooled_pr1.bam pooled_pr2.bam

macs2 callpeak -t pooled_pr1.bam -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n pooled_pr1 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/pooled_pr1_macs2.log
macs2 callpeak -t pooled_pr2.bam -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n pooled_pr2 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/pooled_pr2_macs2.log

sort -k8,8nr peaks/pooled_pr1_peaks.narrowPeak > peaks/pooled_pr1_sorted.narrowPeak
sort -k8,8nr peaks/pooled_pr2_peaks.narrowPeak > peaks/pooled_pr2_sorted.narrowPeak

idr --samples peaks/pooled_pr1_sorted.narrowPeak peaks/pooled_pr2_sorted.narrowPeak \
  --input-file-type narrowPeak --rank p.value \
  --output-file idr/pooled_pseudorep_idr.txt --idr-threshold 0.01 \
  --plot --log-output-file idr/pooled_pseudorep_idr.log

# ─────────────────────────────────────────────────────────────
# STEP 7 — True replicate IDR
# Call peaks on full rep1 and rep2, then run IDR between them.
# This is the most biologically meaningful QC step.
# ─────────────────────────────────────────────────────────────
echo "Step 7: True replicate IDR..."

macs2 callpeak -t $REP1 -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n rep1 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/rep1_macs2.log
macs2 callpeak -t $REP2 -c $CTRL -f BAMPE -g $GENOME_SIZE \
  -n rep2 --outdir peaks/ -p 1e-3 --keep-dup all 2> peaks/rep2_macs2.log

sort -k8,8nr peaks/rep1_peaks.narrowPeak > peaks/rep1_sorted.narrowPeak
sort -k8,8nr peaks/rep2_peaks.narrowPeak > peaks/rep2_sorted.narrowPeak

idr --samples peaks/rep1_sorted.narrowPeak peaks/rep2_sorted.narrowPeak \
  --input-file-type narrowPeak --rank p.value \
  --output-file idr/true_rep_idr.txt --idr-threshold 0.05 \
  --plot --log-output-file idr/true_rep_idr.log

# ─────────────────────────────────────────────────────────────
# STEP 8 — QC ratios
# Rescue ratio < 2 and self-consistency ratio < 2 = PASS
# Optimal peak set: whichever of Np or Nt is larger
# ─────────────────────────────────────────────────────────────
echo "Step 8: Computing QC ratios..."

Nt=$(wc -l < idr/true_rep_idr.txt)
Np=$(wc -l < idr/pooled_pseudorep_idr.txt)
N1=$(wc -l < idr/rep1_pseudorep_idr.txt)
N2=$(wc -l < idr/rep2_pseudorep_idr.txt)

python3 - << PYEOF
Nt=$Nt; Np=$Np; N1=$N1; N2=$N2
rescue = max(Np/Nt, Nt/Np) if Nt > 0 else 0
self_consistency = max(N1/N2, N2/N1) if N2 > 0 else 0
print("")
print("═══════════ QC SUMMARY ═══════════")
print(f"Nt (true rep IDR peaks):      {Nt}")
print(f"Np (pooled pseudorep peaks):  {Np}")
print(f"N1 (rep1 self consistency):   {N1}")
print(f"N2 (rep2 self consistency):   {N2}")
print(f"Rescue ratio:                 {rescue:.3f}  ({'PASS' if rescue < 2 else 'FAIL'})  [threshold < 2]")
print(f"Self consistency ratio:       {self_consistency:.3f}  ({'PASS' if self_consistency < 2 else 'FAIL'})  [threshold < 2]")
print(f"Optimal peak set:             {'Np' if Np > Nt else 'Nt'} ({max(Np,Nt)} peaks)")
print(f"Final peak file:              idr/{'pooled_pseudorep_idr.txt' if Np > Nt else 'true_rep_idr.txt'}")
print("══════════════════════════════════")
PYEOF
