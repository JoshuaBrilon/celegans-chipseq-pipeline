# Genome Setup

## Libuda Lab Users

The N2 genome FASTA and bowtie2 index are already available on Talapas:

```
/projects/libudalab/jbrilon/N2_genome.fasta
/projects/libudalab/jbrilon/N2_genome_index  (.bt2 files)
```

Update `GENOME_INDEX` in `scripts/01_align.sh` to point to this prefix. No further setup needed.

---

## Other Users

### 1. Obtain your genome FASTA

Use your lab's custom *C. elegans* genome FASTA, or download the standard ce11 assembly:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ce11.fa.gz
```

### 2. Build bowtie2 index

Run on a compute node — do not run on the login node:

```bash
srun --partition=interactive --time=120 --mem=8G --pty bash
module load bowtie2/2.4.4
bowtie2-build /path/to/genome.fasta /path/to/genome_index
```

This produces 6 files: `genome_index.1.bt2`, `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, `.rev.2.bt2`

### 3. Update the pipeline

Edit `scripts/01_align.sh` and set:

```bash
GENOME_INDEX="/path/to/your/genome_index"
```

### 4. Control BAM requirement

Your control (input) BAM **must be aligned to the same genome** as your ChIP samples. If using a pre-processed BAM from a public resource:

- Verify chromosome names match: `samtools view -H ctrl.bam | grep "^@SQ"`
- Verify chromosome lengths match those in your ChIP BAMs
- If names differ, use `samtools reheader` to rename chromosomes before running Phase 2

Mismatched control BAMs will cause MACS2 peak calling to fail with "No common chromosome names" error.
