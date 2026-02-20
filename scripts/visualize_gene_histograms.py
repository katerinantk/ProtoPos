#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: visualize_gene_histograms.py.
Description:
    Creates comprehensive visualizations of PRO-seq 3' end distributions
    for the top genes by sequencing depth. Generates:
    1. Individual histogram plots for each gene
    2. Multi-panel comparison plot of all top genes
    3. Summary bar chart of sequencing depths

Dependencies:
    • Python >= 3.10.
    • matplotlib >= 3.5.0.
    • pandas >= 1.3.0.
    • numpy >= 1.20.0.

Usage:
    python scripts/visualize_gene_histograms.py
"""

import os
import sys
import glob

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for WSL.
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
except ImportError:
    print("ERROR: matplotlib required. pip install matplotlib")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas required. pip install pandas")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy required. pip install numpy")
    sys.exit(1)


# Configuration
HISTOGRAM_DIR = "results/histograms"
GENE_DEPTH_FILE = "results/gene_depth.tsv"
OUTPUT_DIR = "results/figures"
TOP_N = 10

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)


def read_histogram_data(histogram_file):
    """
    Read histogram data from text file.

    Args:
        histogram_file (str): Path to histogram file.

    Returns:
        tuple: (metadata_dict, positions_array, distances_array)
    """
    metadata = {}
    positions = []
    distances = []

    with open(histogram_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Parse metadata
            if line.startswith('# Gene:'):
                metadata['gene'] = line.split(':', 1)[1].strip()
            elif line.startswith('# Gene ID:'):
                metadata['gene_id'] = line.split(':', 1)[1].strip()
            elif line.startswith('# Chromosome:'):
                metadata['chrom'] = line.split(':', 1)[1].strip()
            elif line.startswith('# TSS:'):
                metadata['tss'] = int(line.split(':', 1)[1].strip())
            elif line.startswith('# Strand:'):
                metadata['strand'] = line.split(':', 1)[1].strip()
            elif line.startswith('# Gene coordinates:'):
                metadata['coords'] = line.split(':', 1)[1].strip()
            elif line.startswith('#') or line == '':
                continue
            elif line.startswith('position'):
                continue  # Header line
            else:
                # Data line: position, distance_from_TSS
                parts = line.split('\t')
                if len(parts) == 2:
                    positions.append(int(parts[0]))
                    distances.append(int(parts[1]))

    return metadata, np.array(positions), np.array(distances)


def plot_single_gene_histogram(gene_name, histogram_file, gene_depth_info, output_file):
    """
    Create a histogram plot for a single gene.

    Args:
        gene_name (str): Gene symbol.
        histogram_file (str): Path to histogram data file.
        gene_depth_info (dict): Dictionary with gene depth statistics.
        output_file (str): Output PNG file path.
    """
    # Read data
    metadata, positions, distances = read_histogram_data(histogram_file)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), height_ratios=[3, 1])

    # Main histogram: distance from TSS
    ax1.hist(distances, bins=50, color='#2E86AB', alpha=0.7, edgecolor='black')
    ax1.axvline(0, color='red', linestyle='--', linewidth=2, label='TSS', alpha=0.7)
    ax1.set_xlabel('Distance from TSS (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Read Count', fontsize=12, fontweight='bold')
    ax1.set_title(f"{gene_name} - PRO-seq 3' End Distribution",
                  fontsize=14, fontweight='bold', pad=15)
    ax1.legend(fontsize=10)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add statistics text box
    stats_text = (
        f"Total Reads: {len(distances)}\n"
        f"Gene Length: {gene_depth_info['length_kb']:.2f} kb\n"
        f"Depth: {gene_depth_info['depth_per_kb']:.2f} reads/kb\n"
        f"Strand: {metadata['strand']}\n"
        f"Coordinates: {metadata['coords']}"
    )
    ax1.text(0.98, 0.97, stats_text, transform=ax1.transAxes,
             fontsize=9, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Bottom panel: cumulative distribution
    sorted_distances = np.sort(distances)
    cumulative = np.arange(1, len(sorted_distances) + 1) / len(sorted_distances) * 100
    ax2.plot(sorted_distances, cumulative, color='#A23B72', linewidth=2)
    ax2.axvline(0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax2.set_xlabel('Distance from TSS (bp)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cumulative %', fontsize=12, fontweight='bold')
    ax2.set_title('Cumulative Distribution', fontsize=11, fontweight='bold')
    ax2.grid(alpha=0.3, linestyle='--')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def plot_comparison_panel(gene_depth_df, output_file):
    """
    Create a multi-panel comparison plot of top genes.

    Args:
        gene_depth_df (pd.DataFrame): DataFrame with gene depth information.
        output_file (str): Output PNG file path.
    """
    top_genes = gene_depth_df.head(TOP_N)

    # Create figure with grid
    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(4, 3, figure=fig, hspace=0.35, wspace=0.3)

    for idx, (_, row) in enumerate(top_genes.iterrows()):
        gene_name = row['gene_name']
        histogram_file = f"{HISTOGRAM_DIR}/{gene_name}_histogram.txt"

        if not os.path.exists(histogram_file):
            continue

        # Read data
        metadata, positions, distances = read_histogram_data(histogram_file)

        # Determine subplot position
        row_idx = idx // 3
        col_idx = idx % 3
        ax = fig.add_subplot(gs[row_idx, col_idx])

        # Plot histogram
        ax.hist(distances, bins=30, color='#2E86AB', alpha=0.7, edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_xlabel('Distance from TSS (bp)', fontsize=9)
        ax.set_ylabel('Read Count', fontsize=9)
        ax.set_title(f"{gene_name}\n{row['read_count']} reads | "
                     f"{row['depth_per_kb']:.1f} reads/kb",
                     fontsize=10, fontweight='bold')
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Overall title
    fig.suptitle('Top 10 Genes by PRO-seq Sequencing Depth\n3\' End Distribution Patterns',
                 fontsize=16, fontweight='bold', y=0.995)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def plot_depth_comparison(gene_depth_df, output_file):
    """
    Create a bar chart comparing sequencing depths.

    Args:
        gene_depth_df (pd.DataFrame): DataFrame with gene depth information.
        output_file (str): Output PNG file path.
    """
    top_genes = gene_depth_df.head(TOP_N)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Panel 1: Depth per kb
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, TOP_N))
    bars1 = ax1.barh(range(TOP_N), top_genes['depth_per_kb'], color=colors,
                      edgecolor='black', linewidth=1.5)
    ax1.set_yticks(range(TOP_N))
    ax1.set_yticklabels(top_genes['gene_name'], fontsize=11, fontweight='bold')
    ax1.set_xlabel('Sequencing Depth (reads/kb)', fontsize=12, fontweight='bold')
    ax1.set_title('Sequencing Depth per Kilobase', fontsize=13, fontweight='bold', pad=15)
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars1, top_genes['depth_per_kb'])):
        ax1.text(val + 2, i, f'{val:.1f}', va='center', fontsize=9, fontweight='bold')

    # Panel 2: Total read count
    bars2 = ax2.barh(range(TOP_N), top_genes['read_count'], color=colors,
                      edgecolor='black', linewidth=1.5)
    ax2.set_yticks(range(TOP_N))
    ax2.set_yticklabels(top_genes['gene_name'], fontsize=11, fontweight='bold')
    ax2.set_xlabel('Total Read Count', fontsize=12, fontweight='bold')
    ax2.set_title('Total 3\' End Reads per Gene', fontsize=13, fontweight='bold', pad=15)
    ax2.invert_yaxis()
    ax2.grid(axis='x', alpha=0.3, linestyle='--')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars2, top_genes['read_count'])):
        ax2.text(val + 5, i, f'{int(val)}', va='center', fontsize=9, fontweight='bold')

    plt.suptitle('Top 10 Genes - Sequencing Depth Comparison',
                 fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def main():
    """
    Main entry point for visualization.
    """
    print("=" * 80)
    print("PRO-seq Gene Histogram Visualization")
    print("=" * 80)

    # Load gene depth data
    print(f"\nLoading gene depth data from: {GENE_DEPTH_FILE}")
    gene_depth_df = pd.read_csv(GENE_DEPTH_FILE, sep='\t')
    top_genes = gene_depth_df.head(TOP_N)

    print(f"Top {TOP_N} genes:")
    for idx, row in top_genes.iterrows():
        print(f"  {idx+1}. {row['gene_name']}: {row['depth_per_kb']:.2f} reads/kb")

    # Create individual gene plots
    print(f"\n{'='*80}")
    print("Creating individual gene histogram plots...")
    print('='*80)

    for idx, row in top_genes.iterrows():
        gene_name = row['gene_name']
        histogram_file = f"{HISTOGRAM_DIR}/{gene_name}_histogram.txt"

        if not os.path.exists(histogram_file):
            print(f"  WARNING: {histogram_file} not found, skipping")
            continue

        output_file = f"{OUTPUT_DIR}/{gene_name}_histogram.png"
        gene_depth_info = {
            'length_kb': row['length_kb'],
            'depth_per_kb': row['depth_per_kb'],
            'read_count': row['read_count']
        }

        print(f"\n{idx+1}. {gene_name}")
        plot_single_gene_histogram(gene_name, histogram_file, gene_depth_info, output_file)

    # Create comparison panel
    print(f"\n{'='*80}")
    print("Creating multi-gene comparison panel...")
    print('='*80)
    comparison_file = f"{OUTPUT_DIR}/top10_comparison_panel.png"
    plot_comparison_panel(gene_depth_df, comparison_file)

    # Create depth comparison
    print(f"\n{'='*80}")
    print("Creating sequencing depth comparison plot...")
    print('='*80)
    depth_file = f"{OUTPUT_DIR}/sequencing_depth_comparison.png"
    plot_depth_comparison(gene_depth_df, depth_file)

    # Summary
    print(f"\n{'='*80}")
    print("Visualization Complete!")
    print('='*80)
    print(f"\nGenerated files in: {OUTPUT_DIR}/")
    print(f"  - {TOP_N} individual gene histograms")
    print(f"  - 1 multi-panel comparison plot")
    print(f"  - 1 depth comparison chart")
    print(f"\nTotal: {TOP_N + 2} PNG files")
    print('='*80)


if __name__ == "__main__":
    main()
