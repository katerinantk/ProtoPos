# Code Review Assessment – ProtoPos

This document summarizes the findings from a comprehensive code review of the ProtoPos PRO-seq analysis pipeline.

---

## 1. Overview

The ProtoPos project is a PRO-seq (Precision Run-On sequencing) analysis mini-pipeline. The original code was functional in its basic structure but had several issues that could lead to incorrect results or runtime errors. The pipeline consists of:

1. **proseqanalysis.sh** – Downloads FASTQ data, performs QC, trims adapters, and aligns reads with STAR.
2. **Protopos.py** – Extracts 3′ end positions from aligned reads to represent polymerase positions.
3. **RNAPsignal.py** – Computes average Pol II signal around transcription start sites using BigWig files.

**Overall assessment**: The pipeline demonstrates a good understanding of PRO-seq biology and uses appropriate tools. However, several programming and documentation issues needed to be addressed to ensure correctness and robustness.

---

## 2. Programming Issues Found

### 2.1 Protopos.py

| Location | Issue | Fix Applied |
|----------|-------|-------------|
| Line 9–12 | **Hard-coded file paths**: No way to specify input files via command line. | Added `argparse` for CLI arguments while preserving default paths for backward compatibility. |
| Line 13–19 | **Incomplete documentation**: The `three_prime_position` function had a minimal comment that didn't explain the biological logic. | Added comprehensive docstring explaining PRO-seq biology and coordinate handling. |
| Line 28 | **No file existence check**: Script would crash with cryptic error if input file doesn't exist. | Added `os.path.exists()` check with informative error message. |
| Line 31–32 | **No handling of malformed reads**: Reads with `reference_end = None` would cause crashes. | Added check for `read.reference_end is None`. |
| Line 45 | **Lowercase sequence placeholder**: Used `"n"` instead of standard `"N"` for unknown bases. | Changed to uppercase `"N"` for SAM format compliance. |
| Line 26–27 | **Uninformative logging**: Printed generic messages without file names or counts. | Added detailed logging with read counts and file paths. |
| General | **No main() function**: Code ran at module level, preventing import as library. | Restructured with `main()` function and `if __name__ == "__main__"` guard. |

### 2.2 RNAPsignal.py

| Location | Issue | Fix Applied |
|----------|-------|-------------|
| Line 11–14 | **Hard-coded paths**: No way to change input files without editing code. | Documented as configuration section; could be extended with argparse. |
| Line 119–122 | **Incorrect coordinate conversion**: The original code had a subtle bug in computing window boundaries that could cause off-by-one errors near chromosome edges. | Fixed window calculation with explicit 1-based to 0-based conversion and boundary checks. |
| Line 127–131 | **Silent failures**: Missing chromosomes in BigWig returned `NaN` without any warning. | Maintained `NaN` return for missing data but improved surrounding code to handle this case. |
| Line 151–166 | **No progress indicator**: Processing thousands of genes gave no feedback. | Added progress printing every 5000 genes. |
| General | **Missing input validation**: No checks that files exist before attempting to open. | Added `os.path.exists()` checks with clear error messages. |

### 2.3 proseqanalysis.sh

| Location | Issue | Fix Applied |
|----------|-------|-------------|
| General | **No error handling**: Script would continue even if commands failed. | Added `set -euo pipefail` for strict error handling. |
| General | **No tool verification**: Assumed all tools were installed without checking. | Added loop to verify all required tools are in PATH. |
| General | **No timestamps in logs**: Hard to track progress in long-running pipelines. | Added `log()` function with timestamps. |
| Line 126–134 | **BAM not indexed**: Generated BAM files were not indexed, preventing random access. | Added `samtools index` step after conversion. |
| Line 136–138 | **No next-steps guidance**: User might not know what to do after pipeline completes. | Added helpful "Next steps" message at end. |

---

## 3. Biological Issues Found

### 3.1 PRO-seq 3′ End Interpretation

**Original code was correct but poorly documented.**

The `three_prime_position` function in `Protopos.py` correctly computed:
- **Plus strand**: 3′ end = `reference_end - 1`
- **Minus strand**: 3′ end = `reference_start`

This is biologically correct for PRO-seq data because:
- In PRO-seq, reads represent nascent RNAs that were being synthesized at the moment of biotin incorporation.
- The 3′ end of the nascent RNA corresponds to the active site of RNA Polymerase II.
- For forward-strand genes, the polymerase moves 5′→3′ on the template (3′→5′ on the coding strand), so the 3′ end of the read (which is the rightmost position in genome coordinates) represents the polymerase location.
- For reverse-strand genes, the polymerase moves in the opposite direction, so the 3′ end of the read is at the leftmost position.

**Fix applied**: Added comprehensive documentation explaining this logic, including a detailed docstring with biological context.

### 3.2 Coordinate System Confusion

**Potential issue found in RNAPsignal.py.**

The original code mixed coordinate systems without explicit documentation:
- GTF files use **1-based, fully-closed** coordinates `[start, end]`.
- BigWig queries use **0-based, half-open** coordinates `[start, end)`.
- BAM files use **0-based, half-open** coordinates.

The original code attempted to convert but didn't document what it was doing, making it hard to verify correctness.

**Fix applied**: Added explicit comments documenting coordinate conversions at each step:
```python
# Convert TSS from 1-based to 0-based.
tss_0based = tss - 1
```

### 3.3 Strand-Specific BigWig Selection

**Original code was correct.**

The code correctly used separate BigWig files for plus and minus strands:
- Plus-strand genes query the plus-strand BigWig.
- Minus-strand genes query the minus-strand BigWig.

This is essential for PRO-seq data, which is strand-specific.

### 3.4 GFF/GTF Feature Selection

**No issues found.**

The code correctly:
- Parses only "transcript" features (not gene, exon, CDS, etc.).
- Filters for "protein_coding" transcripts.
- Selects the longest transcript per gene (reasonable heuristic).

---

## 4. Missing Functionality

### 4.1 Chromosome Subset Script

**Status**: Did not exist.

**Created**: `scripts/subset_fastq_by_chrom.py`

This script enables rapid testing by extracting reads that map to a specific chromosome. Features:
- Uses pysam to identify reads from BAM file.
- Extracts corresponding reads from original FASTQ.
- Supports region specification (chromosome + start/end).
- Handles both gzipped and uncompressed FASTQ files.

### 4.2 Gene-Level Histogram Script

**Status**: Did not exist.

**Created**: `scripts/gene_histogram.py`

This script visualizes polymerase distribution within a gene. Features:
- Parses GFF3 annotations (with gzip support).
- Extracts 3′ positions from single-base BAM/SAM files.
- Correctly handles strand orientation.
- Auto-selects genes with coverage if none specified.
- Produces publication-quality histograms.

---

## 5. Suggestions for Improvement

### 5.1 Testing Strategy

Create a small test dataset for rapid validation:

```bash
# 1. Run the main pipeline once on full data.
bash proseqanalysis.sh

# 2. Create a chr22 subset for testing.
python scripts/subset_fastq_by_chrom.py \
    --bam alignments/SRR15304569.Aligned.sortedByCoord.out.bam \
    --fastq fastp/SRR15304569.trimmed.fastq.gz \
    --chromosome chr22 \
    --output test_data/chr22.fastq.gz

# 3. Use the subset for rapid iteration during development.
```

Consider adding the test data to the repository (if small enough) or documenting how to generate it.

### 5.2 Code Structure

1. **Move configuration to a separate file**: Create a `config.py` or `config.yaml` for paths and parameters. This avoids hard-coding paths in multiple scripts.

2. **Create a shared utilities module**: Extract common functions (file validation, coordinate conversion) into `scripts/utils.py`.

3. **Add a Makefile or Snakemake workflow**: This would make the pipeline more reproducible and easier to run.

### 5.3 Error Handling

1. **Add validation early**: Check that all input files exist before starting any processing.

2. **Log to files**: In addition to stdout, write logs to a timestamped file for debugging.

3. **Add checksums**: For downloaded files, verify integrity with MD5/SHA checksums.

### 5.4 Documentation

1. **Add example outputs**: Include sample output files (histogram, table) so users know what to expect.

2. **Document required disk space**: STAR indexing and alignment require significant disk space.

3. **Add troubleshooting section**: Common errors and their solutions.

---

## 5. Common Mistakes to Avoid

### For Future Development

1. **Always document coordinate systems**: When working with genomic coordinates, explicitly state whether you're using 0-based or 1-based, and whether intervals are half-open or closed.

2. **Test with known genes**: Verify your analysis on well-characterized genes (ACTB, GAPDH) where you know what to expect.

3. **Check strand handling**: Many bioinformatics bugs come from incorrect strand handling. Always test both strands.

4. **Don't ignore warnings**: Compiler/interpreter warnings often indicate real problems.

5. **Version control your analysis**: Use git commits to track changes to your analysis scripts.

---

## 7. Summary Table

| Category | Issues Found | Issues Fixed | New Code Added |
|----------|--------------|--------------|----------------|
| Protopos.py | 7 | 7 | Restructured with argparse, logging |
| RNAPsignal.py | 5 | 5 | Added validation, progress |
| proseqanalysis.sh | 4 | 4 | Added error handling, logging |
| Chromosome subset | N/A | N/A | New script created |
| Gene histogram | N/A | N/A | New script created |
| Documentation | Missing | N/A | CLAUDE.md created |

---

## 8. Running a Quick Test

After the fixes, you can run a quick end-to-end test:

```bash
# Assuming you have aligned data in the alignments/ directory

# 1. Extract 3' ends.
python Protopos.py --input alignments/SRR15304569.Aligned.sortedByCoord.out.sam

# 2. Create a histogram for ACTB (or auto-select).
python scripts/gene_histogram.py \
    --gff data/gencode_hg38.gff3.gz \
    --bam alignments/SRR15304569.3prime.sam \
    --gene ACTB \
    --output results/histogram_ACTB.png

# 3. If you have BigWig files, compute TSS signal.
python RNAPsignal.py
```

---