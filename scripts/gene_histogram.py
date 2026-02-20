#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: gene_histogram.py.
Description:
    Creates a histogram of 3' end positions within a gene, with TSS at x=0.
    This visualization shows the distribution of RNA Polymerase II positions
    along a selected protein-coding gene.

Dependencies:
    • Python >= 3.10.
    • pysam >= 0.19.0.
    • matplotlib >= 3.5.0.
    • numpy >= 1.20.0.

Usage:
    python scripts/gene_histogram.py --gff annotation.gtf --bam sample.3prime.bam --output hist.png
"""

import argparse
import gzip
import os
import sys

try:
    import pysam
except ImportError:
    print("ERROR: pysam required. pip install pysam")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("ERROR: matplotlib required. pip install matplotlib")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy required. pip install numpy")
    sys.exit(1)


"""Parse GFF3/GTF for protein-coding genes."""


def parse_annotation_for_genes(annot_path):
    """
    Parse a GFF3 or GTF annotation file to extract protein-coding genes.

    This function handles both GFF3 format (key=value attributes) and GTF format
    (key "value" attributes) automatically.

    Args:
        annot_path (str): Path to the annotation file (GFF3 or GTF, optionally gzipped).

    Returns:
        dict: Dictionary mapping gene names to their coordinates and metadata.
              Each value contains: chrom, start, end, strand, gene_id.

    Raises:
        FileNotFoundError: If the annotation file does not exist.
    """
    if not os.path.exists(annot_path):
        raise FileNotFoundError(f"annotation not found: {annot_path}")

    genes = {}

    # Handle gzipped files automatically.
    opener = gzip.open(annot_path, "rt") if annot_path.endswith(".gz") else open(annot_path, "r")

    with opener as f:
        for line in f:
            # Skip comment lines and empty lines.
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields

            # Only process gene features.
            if feature != "gene":
                continue

            # Parse attributes (handle both GFF3 and GTF formats).
            attributes = {}
            for attr in attrs.strip().rstrip(";").split(";"):
                attr = attr.strip()
                if "=" in attr:
                    # GFF3 format: key=value.
                    key, val = attr.split("=", 1)
                    attributes[key] = val
                elif " " in attr:
                    # GTF format: key "value".
                    parts = attr.split(" ", 1)
                    if len(parts) == 2:
                        key = parts[0]
                        val = parts[1].strip('"')
                        attributes[key] = val

            gene_name = attributes.get("gene_name") or attributes.get("Name")
            gene_id = attributes.get("gene_id") or attributes.get("ID")
            gene_type = attributes.get("gene_type") or attributes.get("biotype")

            # Only include protein-coding genes.
            if gene_type != "protein_coding":
                continue

            if not gene_name:
                gene_name = gene_id
            if not gene_name:
                continue

            if gene_name not in genes:
                genes[gene_name] = {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": gene_id,
                }

    return genes


"""Get 3' positions relative to TSS (TSS = 0)."""


def get_3prime_positions_relative_to_tss(bam_path, chrom, start, end, strand):
    """
    Extract 3' end positions from a BAM/SAM file relative to the gene TSS.

    Positions are calculated so that TSS = 0 and positive values indicate
    positions downstream (into the gene body).

    Args:
        bam_path (str): Path to the 3' end BAM/SAM file.
        chrom (str): Chromosome name.
        start (int): Gene start coordinate (1-based, from annotation).
        end (int): Gene end coordinate (1-based, from annotation).
        strand (str): Gene strand ('+' or '-').

    Returns:
        list: List of positions relative to TSS (TSS = 0, positive = downstream).

    Raises:
        FileNotFoundError: If the BAM/SAM file does not exist.
    """
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM/SAM not found: {bam_path}")

    positions = []
    mode = "rb" if bam_path.endswith(".bam") else "r"

    with pysam.AlignmentFile(bam_path, mode) as bam:
        # Return empty list if chromosome not in file.
        if chrom not in bam.references:
            return []

        # For plus strand: TSS is at 'start' (smallest coordinate).
        # For minus strand: TSS is at 'end' (largest coordinate).
        if strand == "+":
            tss = start
        else:
            tss = end

        # Fetch reads in gene region.
        fetch_start = start - 1  # Convert to 0-based.
        fetch_end = end

        try:
            reads = bam.fetch(chrom, fetch_start, fetch_end)
        except ValueError:
            # Fall back to iterating all reads if index not available.
            reads = bam

        for read in reads:
            # Skip unmapped reads.
            if read.is_unmapped:
                continue
            # Skip reads on different chromosome.
            if read.reference_name != chrom:
                continue

            pos = read.reference_start  # 0-based position.

            # Check strand match.
            if strand == "+":
                if read.is_reverse:
                    continue
            else:
                if not read.is_reverse:
                    continue

            # Check if position is within gene bounds.
            if pos < fetch_start or pos >= fetch_end:
                continue

            # Convert to 1-based for calculation.
            pos_1based = pos + 1

            # Calculate position relative to TSS (TSS = 0).
            if strand == "+":
                rel_pos = pos_1based - tss  # Positive = downstream of TSS.
            else:
                rel_pos = tss - pos_1based  # Positive = downstream of TSS (into gene body).

            positions.append(rel_pos)

    return positions


"""Find a gene with sufficient read coverage."""


def find_gene_with_coverage(genes, bam_path, min_reads=10):
    """
    Find a protein-coding gene with at least min_reads mapped reads.

    Tries common housekeeping genes first, then searches all genes.

    Args:
        genes (dict): Dictionary of genes from parse_annotation_for_genes.
        bam_path (str): Path to the 3' end BAM/SAM file.
        min_reads (int): Minimum number of reads required.

    Returns:
        str: Gene name with sufficient coverage, or None if not found.
    """
    # Try common housekeeping genes first.
    priority = ["ACTB", "GAPDH", "MYC", "TP53", "TUBB", "RPL13A"]

    for name in priority:
        if name in genes:
            g = genes[name]
            pos = get_3prime_positions_relative_to_tss(
                bam_path, g["chrom"], g["start"], g["end"], g["strand"]
            )
            if len(pos) >= min_reads:
                return name

    # Try all other genes.
    for name, g in genes.items():
        pos = get_3prime_positions_relative_to_tss(
            bam_path, g["chrom"], g["start"], g["end"], g["strand"]
        )
        if len(pos) >= min_reads:
            return name

    return None


"""Create histogram with TSS at x=0."""


def create_histogram(positions, gene_name, gene_info, n_bins, output_path):
    """
    Create and save a histogram of polymerase positions along a gene.

    The histogram shows the distribution of 3' end positions with TSS at x=0
    and TES marked at the gene length.

    Args:
        positions (list): List of positions relative to TSS.
        gene_name (str): Name of the gene for the title.
        gene_info (dict): Gene metadata (chrom, start, end, strand).
        n_bins (int): Number of histogram bins.
        output_path (str): Path for the output PNG file.
    """
    gene_length = gene_info["end"] - gene_info["start"] + 1

    # Create output directory if needed.
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fig, ax = plt.subplots(figsize=(12, 5))

    # Create histogram if there are positions to plot.
    if len(positions) > 0:
        min_pos = min(positions)
        max_pos = max(positions)
        bins = np.linspace(min_pos, max_pos, n_bins + 1)
        ax.hist(positions, bins=bins, color="steelblue", edgecolor="white", alpha=0.8)

    # Mark TSS and TES.
    ax.axvline(x=0, color="red", linestyle="--", linewidth=2, label="TSS")
    ax.axvline(x=gene_length, color="green", linestyle="--", linewidth=1.5, label="TES")

    ax.set_xlabel("Position relative to TSS (bp)")
    ax.set_ylabel("Number of polymerase positions")
    ax.set_title(
        f"RNA Pol II Distribution: {gene_name} ({gene_info['gene_id']})\n"
        f"{gene_info['chrom']}:{gene_info['start']:,}-{gene_info['end']:,} "
        f"({gene_info['strand']} strand)"
    )

    ax.legend(loc="upper right")

    # Add statistics box.
    stats = f"Reads: {len(positions):,}\nLength: {gene_length:,} bp"
    ax.text(0.02, 0.98, stats, transform=ax.transAxes, verticalalignment="top",
            fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"histogram saved: {output_path}")


"""Main entry point."""


def main():
    """
    Main entry point for the gene histogram script.

    Parses command-line arguments, loads annotation, selects a gene,
    extracts 3' end positions, and creates the histogram.
    """
    parser = argparse.ArgumentParser(description="Plot 3' ends within a gene (TSS at x=0).")
    parser.add_argument("--gff", "-g", required=True, help="GFF3/GTF annotation file.")
    parser.add_argument("--bam", "-b", required=True, help="3' end SAM/BAM file.")
    parser.add_argument("--gene", "-n", default=None, help="Gene name (auto-selects if not given).")
    parser.add_argument("--bins", type=int, default=100, help="Number of bins.")
    parser.add_argument("--output", "-o", required=True, help="Output PNG file.")
    parser.add_argument("--min-reads", type=int, default=10, help="Min reads for auto-selection.")

    args = parser.parse_args()

    print("=" * 60)
    print("Gene Histogram (TSS at x=0)")
    print("=" * 60)

    # Parse annotation.
    print(f"parsing {args.gff}...")
    genes = parse_annotation_for_genes(args.gff)
    print(f"found {len(genes)} protein-coding genes")

    if len(genes) == 0:
        print("ERROR: no genes found")
        sys.exit(1)

    # Select gene.
    if args.gene:
        if args.gene not in genes:
            print(f"ERROR: gene '{args.gene}' not found")
            similar = [g for g in genes if args.gene.upper() in g.upper()][:5]
            if similar:
                print(f"similar: {', '.join(similar)}")
            sys.exit(1)
        gene_name = args.gene
    else:
        print("auto-selecting gene with coverage...")
        gene_name = find_gene_with_coverage(genes, args.bam, args.min_reads)
        if gene_name is None:
            print(f"ERROR: no genes with >= {args.min_reads} reads")
            sys.exit(1)

    gene_info = genes[gene_name]
    print(f"selected: {gene_name}")
    print(f"  {gene_info['chrom']}:{gene_info['start']}-{gene_info['end']} ({gene_info['strand']})")

    # Get positions.
    print("extracting 3' positions...")
    positions = get_3prime_positions_relative_to_tss(
        args.bam,
        gene_info["chrom"],
        gene_info["start"],
        gene_info["end"],
        gene_info["strand"]
    )
    print(f"found {len(positions)} positions")

    # Create histogram.
    print("creating histogram...")
    create_histogram(positions, gene_name, gene_info, args.bins, args.output)

    print("=" * 60)
    print("done")
    print("=" * 60)


if __name__ == "__main__":
    main()
