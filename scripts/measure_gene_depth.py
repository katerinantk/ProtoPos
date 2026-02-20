#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: measure_gene_depth.py.
Description:
    Measures PRO-seq sequencing depth across protein-coding genes.
    For each gene, counts the number of 3' end reads and calculates:
    - Total read count
    - Gene length (in kb)
    - Sequencing depth (reads per kb)

    Outputs a TSV file with genes sorted by sequencing depth.

Dependencies:
    • Python >= 3.10.
    • pysam >= 0.19.0.
    • pandas >= 1.3.0.

Usage:
    python scripts/measure_gene_depth.py \
        --bam alignments/SRR15304570.3prime.bam \
        --gtf data/gencode.v47.annotation.gtf \
        --num-genes 100 \
        --output results/gene_depth.tsv
"""

import argparse
import sys
import os
from collections import defaultdict

try:
    import pysam
except ImportError:
    print("ERROR: pysam is required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas is required. Install with: pip install pandas", file=sys.stderr)
    sys.exit(1)


def parse_gtf_genes(gtf_path, gene_type="protein_coding", num_genes=100):
    """
    Parse GTF file to extract gene coordinates for protein-coding genes.

    Args:
        gtf_path (str): Path to GTF annotation file.
        gene_type (str): Gene biotype to filter (default: "protein_coding").
        num_genes (int): Maximum number of genes to return.

    Returns:
        dict: Dictionary mapping gene_name to (chrom, start, end, strand, gene_id).
    """
    print(f"Parsing GTF file: {gtf_path}")
    print(f"Filtering for gene_type: {gene_type}")

    genes = {}

    with open(gtf_path, 'r') as f:
        for line in f:
            # Skip comments.
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])  # GTF is 1-based.
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Only process 'gene' features.
            if feature != 'gene':
                continue

            # Parse attributes.
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                parts = attr.split(' ', 1)
                if len(parts) == 2:
                    key = parts[0]
                    value = parts[1].strip('"')
                    attr_dict[key] = value

            # Check gene type.
            if 'gene_type' not in attr_dict or attr_dict['gene_type'] != gene_type:
                continue

            # Extract gene name and ID.
            gene_name = attr_dict.get('gene_name', attr_dict.get('gene_id', 'unknown'))
            gene_id = attr_dict.get('gene_id', 'unknown')

            # Store gene info (use gene_name as key for easy access).
            # Convert to 0-based coordinates for pysam.
            genes[gene_name] = {
                'chrom': chrom,
                'start': start - 1,  # Convert to 0-based.
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'length_kb': (end - start + 1) / 1000.0
            }

    print(f"Found {len(genes)} {gene_type} genes in GTF.")
    if num_genes and len(genes) > num_genes:
        # Sort by length descending and keep the top num_genes.
        sorted_names = sorted(genes, key=lambda g: genes[g]['length_kb'], reverse=True)
        genes = {name: genes[name] for name in sorted_names[:num_genes]}
        print(f"Kept top {num_genes} genes by length.")
    return genes


def count_reads_in_genes(bam_path, genes):
    """
    Count PRO-seq 3' end reads within each gene.

    Args:
        bam_path (str): Path to BAM file with 3' end positions.
        genes (dict): Dictionary of gene coordinates from parse_gtf_genes.

    Returns:
        dict: Dictionary mapping gene_name to read count.
    """
    print(f"Counting reads in BAM file: {bam_path}")

    # Check that BAM file exists and is indexed.
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path):
        raise FileNotFoundError(
            f"BAM index not found. Please run: samtools index {bam_path}"
        )

    read_counts = defaultdict(int)

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for gene_name, gene_info in genes.items():
        chrom = gene_info['chrom']
        start = gene_info['start']
        end = gene_info['end']

        # Fetch reads overlapping the gene.
        try:
            for read in bam.fetch(chrom, start, end):
                if not read.is_unmapped:
                    read_counts[gene_name] += 1
        except ValueError:
            # Chromosome not found in BAM (might be chr vs no-chr naming issue).
            # Try alternative naming.
            alt_chrom = chrom.replace('chr', '') if 'chr' in chrom else f'chr{chrom}'
            try:
                for read in bam.fetch(alt_chrom, start, end):
                    if not read.is_unmapped:
                        read_counts[gene_name] += 1
            except ValueError:
                # Still not found, skip this gene.
                continue

    bam.close()

    return read_counts


def calculate_depth_and_sort(genes, read_counts):
    """
    Calculate sequencing depth (reads per kb) and sort genes.

    Args:
        genes (dict): Gene information dictionary.
        read_counts (dict): Read counts per gene.

    Returns:
        pandas.DataFrame: Sorted dataframe with gene depth information.
    """
    data = []

    for gene_name, gene_info in genes.items():
        reads = read_counts.get(gene_name, 0)
        length_kb = gene_info['length_kb']
        depth = reads / length_kb if length_kb > 0 else 0

        data.append({
            'gene_name': gene_name,
            'gene_id': gene_info['gene_id'],
            'chrom': gene_info['chrom'],
            'start': gene_info['start'] + 1,  # Convert back to 1-based for output.
            'end': gene_info['end'],
            'strand': gene_info['strand'],
            'length_kb': round(length_kb, 3),
            'read_count': reads,
            'depth_per_kb': round(depth, 2)
        })

    # Create DataFrame and sort by depth.
    df = pd.DataFrame(data)
    df = df.sort_values('depth_per_kb', ascending=False).reset_index(drop=True)

    return df


def main():
    """
    Main entry point for gene depth measurement.
    """
    parser = argparse.ArgumentParser(
        description="Measure PRO-seq sequencing depth across protein-coding genes."
    )
    parser.add_argument(
        "--bam", "-b",
        required=True,
        help="Path to 3' end BAM file."
    )
    parser.add_argument(
        "--gtf", "-g",
        required=True,
        help="Path to GTF annotation file."
    )
    parser.add_argument(
        "--num-genes", "-n",
        type=int,
        default=100,
        help="Number of genes to analyze (default: 100)."
    )
    parser.add_argument(
        "--gene-type", "-t",
        default="protein_coding",
        help="Gene biotype to filter (default: protein_coding)."
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output TSV file path."
    )

    args = parser.parse_args()

    print("=" * 70)
    print("PRO-seq Gene Depth Analysis")
    print("=" * 70)

    # Step 1: Parse GTF to get gene coordinates.
    try:
        genes = parse_gtf_genes(args.gtf, args.gene_type, args.num_genes)
    except Exception as e:
        print(f"ERROR parsing GTF: {e}", file=sys.stderr)
        sys.exit(1)

    if len(genes) == 0:
        print("ERROR: No genes found matching criteria.", file=sys.stderr)
        sys.exit(1)

    # Step 2: Count reads in each gene.
    print()
    try:
        read_counts = count_reads_in_genes(args.bam, genes)
    except Exception as e:
        print(f"ERROR counting reads: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Calculate depth and sort.
    print("\nCalculating sequencing depth...")
    df = calculate_depth_and_sort(genes, read_counts)

    # Step 4: Write output.
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df.to_csv(args.output, sep='\t', index=False)
    print(f"\nResults saved to: {args.output}")

    # Print summary.
    print("\n" + "=" * 70)
    print("Top 10 Genes by Sequencing Depth")
    print("=" * 70)
    print(df.head(10).to_string(index=False))

    print("\n" + "=" * 70)
    print("Summary Statistics")
    print("=" * 70)
    print(f"Total genes analyzed:     {len(df)}")
    print(f"Genes with reads:         {(df['read_count'] > 0).sum()}")
    print(f"Mean depth (reads/kb):    {df['depth_per_kb'].mean():.2f}")
    print(f"Median depth (reads/kb):  {df['depth_per_kb'].median():.2f}")
    print(f"Max depth (reads/kb):     {df['depth_per_kb'].max():.2f}")
    print("=" * 70)


if __name__ == "__main__":
    main()
