# CLAUDE.md – ProtoPos

This file provides **mandatory rules** and **project context** for Claude Code when editing or generating code in ProtoPos.

---

## **1. Project Overview**

**ProtoPos** is a nascent-RNA sequencing analysis mini-pipeline that extracts RNA Polymerase II active-site positions from aligned BAM files. It supports both **PRO-seq** (antisense library, 5' end of read = active site) and **GRO-seq** (sense library, 3' end of read = active site). The pipeline performs:

- **Read acquisition and QC**: Downloads PRO-seq FASTQ files from SRA, runs FastQC, trims adapters with fastp.
- **UMI handling**: Extracts UMIs with `umi_tools extract` and removes PCR duplicates with `umi_tools dedup`.
- **Genome alignment**: Maps trimmed reads with bowtie2 (splice-unaware, appropriate for nascent unspliced RNA).
- **Pol II position extraction**: Collapses each read to a single base at the RNA Polymerase II active site. For PRO-seq: 5' end of R1 = active site. For GRO-seq: 3' end of the sense read = active site.
- **Signal analysis**: Computes and visualises RNAPII signal around transcription start sites using BigWig coverage files.

**ProtoPos is a standalone analysis tool** within the broader transcription modelling framework:

```
../POLYARIS/   ← Transcription simulation engine
../KINEMA/     ← Analytical Markov chain kinetics
../GENEFIT/    ← CMA-ES parameter fitting
../SCAN/       ← Statistical analyses and figures
../ProtoPos/   ← You are here (PRO-seq data processing)
```

---

## **2. Directory Structure**

```
ProtoPos/
├── Protopos.py              # Extract Pol II positions from aligned BAMs
├── RNAPsignal.py            # Compute RNAPII signal around TSS from BigWig
├── proseqanalysis.sh        # Main bash pipeline: download, QC, trim, align, dedup
├── scripts/                 # Utility scripts
├── data/                    # Input data (genome FASTA, GTF annotation, BigWig)
├── fastq/                   # Raw FASTQ files (created by pipeline)
├── umi/                     # UMI-extracted FASTQ files (created by pipeline)
├── fastp/                   # Trimmed FASTQ files (created by pipeline)
├── qc/                      # Quality control reports (FastQC, fastp, dedup stats)
├── alignments/              # Bowtie2 BAMs: sorted, deduplicated, 3'-end
├── bowtie2_index/           # Bowtie2 genome index
├── results/                 # Analysis outputs: figures, tables
├── .gitignore
├── README.md
└── CLAUDE.md                # This file
```

---

## **3. How to Run**

### **3.1 Full Pipeline (Single-End PRO-seq)**

```bash
# Activate the conda environment.
conda activate protopos

# Place reference files in data/.
#   data/genome.fa                    – Reference genome FASTA.
#   data/gencode.v47.annotation.gtf   – Gene annotation GTF.

# Run the main pipeline.
bash proseqanalysis.sh
```

This downloads FASTQs, extracts UMIs, trims adapters, aligns with bowtie2, and deduplicates.

### **3.2 Extract Pol II Positions**

```bash
python Protopos.py
```

Reads deduplicated BAMs from `alignments/` and writes `*.3prime.sorted.bam` files where each read is a single base at the Pol II active site.

### **3.3 Compute TSS Signal from BigWig**

```bash
python RNAPsignal.py
```

Reads annotation and BigWig files, computes mean RNAPII signal in a ±200 bp window around each TSS, and produces a violin plot.

---

## **4. Environment Expectations**

Claude must always assume:

- **You are running inside WSL** (Linux environment) or Windows.
- The **conda environment `protopos`** is available and contains: bowtie2, samtools, fastp, fastqc, sra-tools, umi_tools, bedtools, pysam, pyBigWig, pandas, numpy, matplotlib.
- Activate with: `source ~/miniconda3/bin/activate protopos`

---

## **5. Instructions for Claude**

### **5.1 General Behavior**

- Prefer **minimal, local, safe edits** that preserve existing structure.
- **Do not** attempt large-scale rewrites or architectural changes unless explicitly asked.
- `Protopos.py` is both a **standalone script** and an **importable module**. The SCAN project imports `pol2_position` and `extract_active_sites` from it.

### **5.2 Coding Style**

#### **Comments & Documentation**
- All explanatory comments must be full sentences ending with a **full stop**.
- Function-level docstrings must be **Google-style**, including:

```python
"""Short summary.

Args:
    param_name (type): Description.
Returns:
    type: Description.
Raises:
    ErrorType: Description.
"""
```

- Include parameter types, return types, assumptions, and biological meaning.

#### **Titles / Subtitles**
- Titles must use: `"""Title"""`
- Subtitles must use: `'''Subtitle'''`
- No alternative formats.

#### **Code Quality**
- Prefer **small, testable functions** rather than long monolithic blocks.
- Strive for optimized efficient code that matches the style of the rest of the file.

### **5.3 Biological Correctness**

This is the most critical aspect of ProtoPos. All code must be biologically correct:

- **Pol II position (PRO-seq)**: The **5' end of the sequencing read** = 3' end of nascent RNA = Pol II active site. Forward alignment: `reference_start`. Reverse alignment: `reference_end - 1`.
- **Pol II position (GRO-seq)**: The **3' end of the sense read** = Pol II active site. Forward alignment: `reference_end - 1`. Reverse alignment: `reference_start`.
- **Strandedness**: Both PRO-seq and GRO-seq are strand-specific. Always respect strand when extracting or analysing signal.
- **Coordinate conventions**: BAM uses 0-based coordinates; GTF uses 1-based. Be explicit about conversions.
- **UMI handling**: Some PRO-seq protocols include UMIs (e.g. Mahat et al. 2016). Check read structure (base composition at position 1) to confirm. If position 1 shows strong G-bias (~60%), there are no 5' UMIs.
- **Aligner choice**: PRO-seq reads nascent **unspliced** RNA. Use splice-unaware aligners (bowtie2), not splice-aware ones (STAR, HISAT2).
- Never alter biological or bioinformatic behaviour silently.

### **5.4 Error Handling**

- Use informative `print()` messages with sample IDs and file paths.
- Check for missing input files before processing.
- In bash scripts, validate that expected output files exist before proceeding.

### **5.5 Testing**

- For quick tests, use **chr22** as the reference (~50 Mb, many genes, fast).
- Use 500K–1M reads for rapid iteration.
- Verify Pol II extraction by checking that forward reads use `reference_start` and reverse reads use `reference_end - 1`.

---

## **6. Author Header Block**

For every **Python file** that Claude creates or substantially modifies, include:

```python
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: <filename>.py.
Description:
    <Brief description of what this file does.>

Dependencies:
    • Python >= 3.10, Python <= 3.13.
    • <list relevant dependencies>

Usage:
    <how to run this script>
"""
```

---

## **7. Git Policy**

- Never git add, commit, stash, tag, or pull.
- Only update and deal with the local version of this tool.
- Never leave comments related to version changes like: `# Version: 1.1 — removed X`.
- Use `.gitkeep` files to preserve empty directories in git.

---

## **8. Intermediate Files**

- **Raw FASTQ**: `fastq/<SRA_ID>.fastq.gz`
- **UMI-extracted FASTQ**: `umi/<SRA_ID>.umi.fastq.gz`
- **Trimmed FASTQ**: `fastp/<SRA_ID>.trimmed.fastq.gz`
- **QC reports**: `qc/<SRA_ID>.fastp.html`, `qc/<SRA_ID>.fastp.json`
- **Aligned BAM**: `alignments/<SRA_ID>.sorted.bam`
- **Deduplicated BAM**: `alignments/<SRA_ID>.dedup.bam`
- **Pol II BAM**: `alignments/<SRA_ID>.3prime.sorted.bam`

---

## **9. Tips & Tricks**

### **Working with This Codebase**

- **Disk space**: BAM and FASTQ files can be large. Use chr22 subsets for development.
- **RAM**: Bowtie2 uses ~3.2 GB for the full human genome index, ~200 MB for chr22 alone.
- **UMI detection**: Check base composition at position 1 of raw reads. Uniform = UMI present; strong G-bias = no UMI (biotin run-on signal).
- **Adapter**: The standard PRO-seq 3' adapter is `TGGAATTCTCGGGTGCCAAGG` (Illumina TruSeq small RNA). Confirm by checking fastp adapter detection statistics.

### **Working with Claude Code Effectively**

- Use subagents for broad file exploration.
- Background long downloads or alignment runs with `run_in_background`.
- Always activate the conda environment before running tools: `source ~/miniconda3/bin/activate protopos`.

---

## **10. End of Instructions**

Claude must follow all rules above.
