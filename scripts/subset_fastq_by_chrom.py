#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: subset_fastq_by_chrom.py.
Description:
    Subsets FASTQ files to only include reads that map to a specific chromosome
    or genomic region. This is useful for creating small test datasets for rapid
    pipeline validation.

    The script works by:
    1. Reading a BAM file to identify which reads map to the target chromosome.
    2. Extracting the corresponding reads from the original FASTQ file(s).
    3. Writing the subset to new FASTQ file(s).

    This approach requires that alignment has already been performed, but ensures
    that only reads with valid mappings are included.

Dependencies:
    • Python >= 3.10.
    • pysam >= 0.19.0.

Usage:
    # Single-end data:
    python scripts/subset_fastq_by_chrom.py \\
        --bam alignments/sample.sorted.bam \\
        --fastq fastp/sample.trimmed.fastq.gz \\
        --chromosome chr22 \\
        --output test_chr22.fastq.gz

    # Optionally specify a region:
    python scripts/subset_fastq_by_chrom.py \\
        --bam alignments/sample.sorted.bam \\
        --fastq fastp/sample.fastq.gz \\
        --chromosome chr22 \\
        --start 20000000 \\
        --end 30000000 \\
        --output test_chr22_region.fastq.gz
"""

import argparse
import gzip
import os
import sys
from typing import Set, Optional

try:
    import pysam
except ImportError:
    print("ERROR: pysam is required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)


def get_read_names_from_bam(
    bam_path: str,
    chromosome: str,
    start: Optional[int] = None,
    end: Optional[int] = None
) -> Set[str]:
    """
    Extract read names from a BAM file for reads mapping to a specific region.

    This function opens a BAM file and extracts all read names (query names) for
    reads that map to the specified chromosome, or optionally a specific region
    within that chromosome.

    Args:
        bam_path (str): Path to the BAM file (must be indexed with .bai file).
        chromosome (str): Chromosome name (e.g., "chr22").
        start (int, optional): Start position of region (1-based, inclusive).
        end (int, optional): End position of region (1-based, inclusive).

    Returns:
        Set[str]: Set of read names mapping to the specified region.

    Raises:
        FileNotFoundError: If BAM file or index not found.
        ValueError: If chromosome not found in BAM header.
    """
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    # Check for BAM index.
    bai_path = bam_path + ".bai"
    csi_path = bam_path + ".csi"
    if not os.path.exists(bai_path) and not os.path.exists(csi_path):
        raise FileNotFoundError(
            f"BAM index not found. Please run: samtools index {bam_path}"
        )

    read_names = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # Verify chromosome exists in BAM header.
        if chromosome not in bam.references:
            available = ", ".join(list(bam.references)[:10])
            raise ValueError(
                f"Chromosome '{chromosome}' not found in BAM. "
                f"Available chromosomes include: {available}..."
            )

        # Fetch reads from the specified region.
        # pysam.fetch uses 0-based coordinates, but we accept 1-based from user.
        if start is not None and end is not None:
            # Convert 1-based to 0-based for pysam.
            region_start = start - 1
            region_end = end
            print(f"Fetching reads from {chromosome}:{start}-{end}...")
        else:
            # Fetch entire chromosome.
            region_start = None
            region_end = None
            print(f"Fetching reads from entire chromosome {chromosome}...")

        for read in bam.fetch(chromosome, region_start, region_end):
            # Skip unmapped reads (though they shouldn't appear in a fetch).
            if read.is_unmapped:
                continue
            read_names.add(read.query_name)

    print(f"Found {len(read_names)} unique read names mapping to {chromosome}.")
    return read_names


def open_fastq(path: str, mode: str = "rt"):
    """
    Open a FASTQ file, handling gzip compression automatically.

    Args:
        path (str): Path to the FASTQ file.
        mode (str): File open mode ('rt' for read text, 'wt' for write text).

    Returns:
        file handle: Open file handle for reading or writing.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    else:
        return open(path, mode)


def subset_fastq(
    fastq_in: str,
    fastq_out: str,
    read_names: Set[str]
) -> tuple:
    """
    Subset a FASTQ file to only include reads with names in the given set.

    FASTQ format:
    - Line 1: @read_name (header line starting with @)
    - Line 2: Sequence
    - Line 3: + (optionally followed by read name again)
    - Line 4: Quality scores

    Args:
        fastq_in (str): Path to input FASTQ file.
        fastq_out (str): Path to output FASTQ file.
        read_names (Set[str]): Set of read names to keep.

    Returns:
        tuple: (total_reads, kept_reads) counts.
    """
    if not os.path.exists(fastq_in):
        raise FileNotFoundError(f"Input FASTQ not found: {fastq_in}")

    # Create output directory if needed.
    out_dir = os.path.dirname(fastq_out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    total_reads = 0
    kept_reads = 0

    # Determine if output should be compressed.
    write_mode = "wt"

    with open_fastq(fastq_in, "rt") as fin, open_fastq(fastq_out, write_mode) as fout:
        while True:
            # Read all four lines of a FASTQ record.
            header = fin.readline()
            if not header:
                break  # End of file.

            sequence = fin.readline()
            plus_line = fin.readline()
            quality = fin.readline()

            # Validate FASTQ format.
            if not header.startswith("@"):
                print(f"WARNING: Malformed FASTQ header: {header.strip()}", file=sys.stderr)
                continue

            total_reads += 1

            # Extract read name from header.
            # Header format: @read_name optional_description
            # The read name is everything after @ up to the first whitespace.
            read_name = header[1:].split()[0]

            # Check if this read should be kept.
            if read_name in read_names:
                fout.write(header)
                fout.write(sequence)
                fout.write(plus_line)
                fout.write(quality)
                kept_reads += 1

            # Progress indicator.
            if total_reads % 1000000 == 0:
                print(f"  Processed {total_reads:,} reads, kept {kept_reads:,}...")

    return (total_reads, kept_reads)


def main():
    """
    Main entry point for chromosome-based FASTQ subsetting.
    """
    parser = argparse.ArgumentParser(
        description="Subset FASTQ files to reads mapping to a specific chromosome.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Subset to all reads mapping to chr22:
    python subset_fastq_by_chrom.py \\
        --bam sample.bam --fastq sample.fastq.gz \\
        --chromosome chr22 --output chr22_test.fastq.gz

    # Subset to a specific region:
    python subset_fastq_by_chrom.py \\
        --bam sample.bam --fastq sample.fastq.gz \\
        --chromosome chr22 --start 20000000 --end 30000000 \\
        --output chr22_region.fastq.gz
        """
    )

    parser.add_argument(
        "--bam", "-b",
        required=True,
        help="Path to indexed BAM file."
    )
    parser.add_argument(
        "--fastq", "-f",
        required=True,
        help="Path to input FASTQ file (may be gzipped)."
    )
    parser.add_argument(
        "--chromosome", "-c",
        required=True,
        help="Target chromosome (e.g., 'chr22')."
    )
    parser.add_argument(
        "--start", "-s",
        type=int,
        default=None,
        help="Start position of region (1-based, optional)."
    )
    parser.add_argument(
        "--end", "-e",
        type=int,
        default=None,
        help="End position of region (1-based, optional)."
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path to output FASTQ file."
    )

    args = parser.parse_args()

    # Validate region arguments.
    if (args.start is None) != (args.end is None):
        print("ERROR: Both --start and --end must be specified together, or neither.",
              file=sys.stderr)
        sys.exit(1)

    if args.start is not None and args.start >= args.end:
        print("ERROR: --start must be less than --end.", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("FASTQ Chromosome Subset")
    print("=" * 60)
    print(f"BAM file:    {args.bam}")
    print(f"FASTQ file:  {args.fastq}")
    print(f"Chromosome:  {args.chromosome}")
    if args.start is not None:
        print(f"Region:      {args.start:,} - {args.end:,}")
    print(f"Output:      {args.output}")
    print()

    # Step 1: Get read names from BAM.
    print("Step 1: Extracting read names from BAM...")
    try:
        read_names = get_read_names_from_bam(
            args.bam,
            args.chromosome,
            args.start,
            args.end
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if len(read_names) == 0:
        print("WARNING: No reads found mapping to the specified region.", file=sys.stderr)
        print("Output file will be empty.")

    # Step 2: Subset FASTQ.
    print("\nStep 2: Subsetting FASTQ file...")
    try:
        total, kept = subset_fastq(args.fastq, args.output, read_names)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Report results.
    print()
    print("=" * 60)
    print("Results")
    print("=" * 60)
    print(f"Total reads in FASTQ: {total:,}")
    print(f"Reads kept:           {kept:,}")
    if total > 0:
        print(f"Percentage kept:      {100 * kept / total:.2f}%")
    print(f"Output file:          {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
