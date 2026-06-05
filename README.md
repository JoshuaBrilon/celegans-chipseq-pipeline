# C. elegans ChIP-seq Pipeline

A two-phase paired-end ChIP-seq alignment and peak calling pipeline for *C. elegans*, developed in the Libuda Lab at the University of Oregon. Validated using ENCODE experiment [ENCSR512EIF](https://www.encodeproject.org/experiments/ENCSR512EIF/) — HSF-1::eGFP ChIP-seq in *C. elegans*.

---

## Why You Need a Control BAM


A ChIP-seq experiment pulls down DNA bound by your protein of interest using an antibody. But some DNA regions are naturally more accessible or "sticky" and get pulled down regardless of your protein. Without a control, you cannot tell whether a peak is real binding or just background noise.

The **control (input) sample** is prepared identically to your ChIP sample but **without the antibody pulldown step**. It captures the background chromatin landscape of your organism. MACS2 uses this to subtract background and only call peaks where your protein is genuinely enriched above what you would expect by chance.


The control must be:
- Sequenced from the same experiment as your ChIP samples
- Aligned to the **exact same genome** as your ChIP samples
- Processed through the same alignment pipeline (`00_create_ctrl_bam.sh`)

See `scripts/00_create_ctrl_bam.sh` for full instructions.

---

## Overview

**Script 00** — Process your control/input FASTQs into a control BAM  
**Phase 1** — Process each ChIP replicate: trim → align → deduplicate → filter  
**Phase 2** — Signal generation, peak calling, pseudoreplication, IDR QC

### Pipeline Diagram

```
Input FASTQs (ctrl)          ChIP FASTQs (rep1, rep2)
       │                           │
       ▼                           ▼
00_create_ctrl_bam.sh         01_align.sh
       │                           │
       ▼                           ▼
ctrl_filtered.bam      rep1_filtered.bam / rep2_filtered.bam
       │                           │
       └──────────┬────────────────┘
                  ▼
            02_phase2.sh
                  │
       ┌──────────┼──────────────┐
       ▼          ▼              ▼
  signal/      peaks/          idr/
  bigWigs    narrowPeaks    IDR peaks + QC
```

---

## Requirements

### System
- Linux (tested on SLURM HPC — University of Oregon Talapas cluster)
- SLURM job scheduler

### Tools

**System modules (loaded via `module load`):**
- `samtools/1.19`
- `bowtie2/2.4.4`
- `cutadapt/4.9`

**Conda environment setup:**
```bash
conda create -n chipseq python=3.10 -y
conda activate chipseq
conda install -c bioconda -c conda-forge picard macs2 deeptools sambamba idr -y
pip install "numpy<1.24"   # required for idr compatibility
```

---

## Genome Setup

> **Libuda Lab users:** The N2 genome FASTA and bowtie2 index are already available on Talapas:
> ```
> /projects/libudalab/jbrilon/N2_genome.fasta
> /projects/libudalab/jbrilon/N2_genome_index
> ```
> Update `GENOME_INDEX` in `scripts/01_align.sh` and `scripts/00_create_ctrl_bam.sh` to point to this prefix. No further genome setup needed.

For other users, see `docs/genome_setup.md` for instructions on building your own index.

---

## Input Files

| File | Description | How to create |
|------|-------------|---------------|
| `ctrl_PE1.fastq.gz` | Control/input, read 1 | From sequencing core |
| `ctrl_PE2.fastq.gz` | Control/input, read 2 | From sequencing core |
| `rep1_PE1.fastq.gz` | Replicate 1 ChIP, read 1 | From sequencing core |
| `rep1_PE2.fastq.gz` | Replicate 1 ChIP, read 2 | From sequencing core |
| `rep2_PE1.fastq.gz` | Replicate 2 ChIP, read 1 | From sequencing core |
| `rep2_PE2.fastq.gz` | Replicate 2 ChIP, read 2 | From sequencing core |

---

## Running the Pipeline

### Step 0 — Create control BAM (required)

**Do this before anything else.** Your control FASTQs must be processed through the same alignment pipeline as your ChIP samples so everything is aligned to the same genome.

```bash
sbatch slurm/submit_ctrl.sh ctrl_PE1.fastq.gz ctrl_PE2.fastq.gz
```

This produces `ctrl_filtered.bam`. See `scripts/00_create_ctrl_bam.sh` for a full explanation of what the control is and why it is required.

### Step 1 — Align ChIP replicates

```bash
sbatch slurm/submit_align.sh rep1_PE1.fastq.gz rep1_PE2.fastq.gz rep1
sbatch slurm/submit_align.sh rep2_PE1.fastq.gz rep2_PE2.fastq.gz rep2
```

These can be submitted at the same time and will run in parallel. Each produces a `{sample}_filtered.bam`.

### Step 2 — Run Phase 2

Once `ctrl_filtered.bam`, `rep1_filtered.bam`, and `rep2_filtered.bam` are all ready:

```bash
sbatch slurm/submit_phase2.sh
```

**Outputs:**
```
signal/
  rep1_foldchange.bigwig     ← ChIP vs control fold change (rep1)
  rep1_pvalue.bigwig         ← p-value track (rep1)
  rep2_foldchange.bigwig
  rep2_pvalue.bigwig
  pooled_foldchange.bigwig   ← pooled rep1+rep2 vs control
  pooled_pvalue.bigwig
peaks/
  rep1_peaks.narrowPeak      ← MACS2 peaks per replicate
  rep2_peaks.narrowPeak
  rep1_pr{1,2}_peaks.narrowPeak  ← pseudoreplicate peaks
  rep2_pr{1,2}_peaks.narrowPeak
  pooled_pr{1,2}_peaks.narrowPeak
idr/
  true_rep_idr.txt           ← Nt: peaks reproducible between rep1 and rep2
  pooled_pseudorep_idr.txt   ← Np: peaks from pooled pseudoreplicates
  rep1_pseudorep_idr.txt     ← N1: rep1 self-consistency
  rep2_pseudorep_idr.txt     ← N2: rep2 self-consistency
  *.png                      ← IDR diagnostic plots
```

---

## QC Interpretation

Phase 2 prints a QC summary at the end:

| Metric | Formula | Pass threshold | Meaning |
|--------|---------|----------------|---------|
| Rescue ratio | max(Np/Nt, Nt/Np) | < 2 | Are pooled pseudoreps consistent with true reps? |
| Self-consistency ratio | max(N1/N2, N2/N1) | < 2 | Are the two replicates consistent with each other? |

**Optimal peak set:** whichever is larger between Np and Nt. This is your final list of binding sites for downstream analysis.

**Example output (ENCODE HSF-1 validation run):**
```
═══════════ QC SUMMARY ═══════════
Nt (true rep IDR peaks):      743
Np (pooled pseudorep peaks):  799
N1 (rep1 self consistency):   516
N2 (rep2 self consistency):   634
Rescue ratio:                 1.075  (PASS)  [threshold < 2]
Self consistency ratio:       1.229  (PASS)  [threshold < 2]
Optimal peak set:             Np (799 peaks)
Final peak file:              idr/pooled_pseudorep_idr.txt
══════════════════════════════════
```

---

## Visualization in IGV

Load results in [IGV Desktop](https://igv.org/doc/desktop/):

1. **Genomes** → **Load Genome from File** → select your genome FASTA
2. **File** → **Load from File** → select `idr/pooled_pseudorep_idr.txt` (optimal peak set)
3. **File** → **Load from File** → select `signal/rep1_foldchange.bigwig` (fold change signal track)
4. Navigate to a peak region to visualize binding

> **Note for ENCODE test data:** The example IGV screenshots in `example_output/` use the ENCODE HSF-1 dataset. The top peak is at `II_1:11784246-11784803`.

---

## ENCODE Test Data

This pipeline was validated using ENCODE experiment [ENCSR512EIF](https://www.encodeproject.org/experiments/ENCSR512EIF/).

```bash
# Replicate 1
wget "https://www.encodeproject.org/files/ENCFF942RFW/@@download/ENCFF942RFW.fastq.gz" -O rep1_PE1.fastq.gz
wget "https://www.encodeproject.org/files/ENCFF526RUB/@@download/ENCFF526RUB.fastq.gz" -O rep1_PE2.fastq.gz

# Replicate 2
wget "https://www.encodeproject.org/files/ENCFF629MDQ/@@download/ENCFF629MDQ.fastq.gz" -O rep2_PE1.fastq.gz
wget "https://www.encodeproject.org/files/ENCFF046ECO/@@download/ENCFF046ECO.fastq.gz" -O rep2_PE2.fastq.gz
```

> **Note on ENCODE control:** ENCFF089TNG (the ENCODE control BAM for this experiment) was aligned to the standard ce11 genome and has different chromosome names and lengths than the Libuda Lab custom N2 genome. For the validation run, chromosome names were remapped using `samtools reheader`. For your own experiments, always create your control BAM using `00_create_ctrl_bam.sh` with FASTQs aligned to your own genome — this is the correct approach and avoids all mismatch issues.

---

## Pipeline Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Adapter | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | Illumina TruSeq |
| Min base quality | Q20 | cutadapt |
| Min read length | 20bp | cutadapt |
| MAPQ threshold | ≥30 | sambamba filter |
| Proper pairs only | yes | sambamba filter |
| Remove mitochondrial | yes | sambamba filter |
| Genome size | 9e7 | *C. elegans* effective genome size for MACS2 |
| MACS2 p-value | 1e-3 | loose threshold — IDR does the final filtering |
| IDR threshold (pseudorep) | 0.05 | |
| IDR threshold (pooled pseudorep) | 0.01 | |
| IDR threshold (true rep) | 0.05 | |

---

## Repository Structure

```
celegans-chipseq-pipeline/
├── README.md
├── .gitignore
├── scripts/
│   ├── 00_create_ctrl_bam.sh   ← process control FASTQs (run first)
│   ├── 01_align.sh             ← process ChIP FASTQs
│   └── 02_phase2.sh            ← peak calling and IDR
├── slurm/
│   ├── submit_ctrl.sh          ← SLURM wrapper for script 00
│   ├── submit_align.sh         ← SLURM wrapper for script 01
│   └── submit_phase2.sh        ← SLURM wrapper for script 02
├── docs/
│   └── genome_setup.md         ← genome index instructions
└── example_output/
    └── (IGV screenshots)
```

---

## Citations

- [ENCODE ChIP-seq guidelines](https://www.encodeproject.org/chip-seq/transcription_factor/)
- [IDR framework](https://github.com/nboley/idr)
- [MACS2](https://github.com/macs3-project/MACS)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/)
- [deepTools](https://deeptools.readthedocs.io/)
- [Picard](https://broadinstitute.github.io/picard/)
- [sambamba](https://lomereiter.github.io/sambamba/)

---

## Contact

Developed by Joshua Brilon — Libuda Lab, University of Oregon  
For questions about the Libuda Lab genome/index, contact the lab directly.
