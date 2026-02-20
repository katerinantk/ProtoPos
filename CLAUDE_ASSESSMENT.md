# Code Review Assessment – ProtoPos

A review of the ProtoPos PRO-seq analysis pipeline.

---

## 1. Overview

ProtoPos is a PRO-seq mini-pipeline that takes FASTQ files, aligns them to a genome, and extracts single-nucleotide 3′ end positions representing RNA polymerase locations. **The core code is correct and demonstrates solid understanding of PRO-seq biology.**

---

## 2. What the Student Did Well

### Correct 3′ End Logic
The `three_prime_position()` function correctly handles strand orientation:
- **Plus strand**: `reference_end - 1` (rightmost base)
- **Minus strand**: `reference_start` (leftmost base)

This is exactly right—the student clearly understands that the 3′ end of nascent RNA marks the polymerase active site.

### Clean SAM Output
Creating new single-base alignments with `cigarstring = "1M"` is the right approach. Important fields (flags, mapping quality, reference ID) are properly preserved.

### Good Pipeline Structure
The bash pipeline follows a sensible workflow: Download → QC → Trim → Align → Convert. Tool choices (STAR, fastp) are appropriate.

### Appropriate Use of pysam
Correct use of `AlignmentFile` for reading/writing, with proper header transfer to output.

---

## 3. Minor Issues Fixed

| File | Issue | Notes |
|------|-------|-------|
| `Protopos.py` | Used lowercase `"n"` for placeholder sequence | Changed to `"N"` (SAM convention) |
| `proseqanalysis.sh` | No error handling | Added `set -euo pipefail` so failures don't go unnoticed |

---

## 4. Enhancements Added

These weren't bugs—just convenience improvements:

- Added `argparse` to `Protopos.py` for flexibility (defaults still work)
- Added input file checks with clear error messages
- Added timestamps and progress logging to bash script
- Added `samtools index` step for BAM files
- Created helper scripts for testing (`subset_fastq_by_chrom.py`, `gene_histogram.py`)

---

## 5. Quick Test

```bash
# Extract 3' ends:
python Protopos.py

# Or with explicit input:
python Protopos.py --input alignments/SRR15304569.Aligned.sortedByCoord.out.sam
```

---

*Good work—the biology and core implementation are correct.*
