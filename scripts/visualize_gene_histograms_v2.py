#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: visualize_gene_histograms_v2.py.
Description:
    Improved visualizations for sparse PRO-seq data.
    Changes from v1:
    - Removed useless cumulative distribution panel
    - Larger bins to reduce noise
    - Smoothed density overlay
    - Better handling of low-coverage data

Dependencies:
    • Python >= 3.10.
    • matplotlib >= 3.5.0.
    • pandas >= 1.3.0.
    • numpy >= 1.20.0.
    • scipy >= 1.7.0.

Usage:
    python scripts/visualize_gene_histograms_v2.py
"""

import os
import sys

try:
    import matplotlib
    matplotlib.use('Agg')
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

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.interpolate import interp1d
    SCIPY_AVAILABLE = True
except ImportError:
    print("WARNING: scipy not available, will skip smoothing")
    SCIPY_AVAILABLE = False


# Configuration
HISTOGRAM_DIR = "results/histograms"
GENE_DEPTH_FILE = "results/gene_depth.tsv"
OUTPUT_DIR = "results/figures"
TOP_N = 10

os.makedirs(OUTPUT_DIR, exist_ok=True)


def read_histogram_data(histogram_file):
    """Read histogram data from text file."""
    metadata = {}
    positions = []
    distances = []

    with open(histogram_file, 'r') as f:
        for line in f:
            line = line.strip()
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
            elif line.startswith('#') or line == '' or line.startswith('position'):
                continue
            else:
                parts = line.split('\t')
                if len(parts) == 2:
                    positions.append(int(parts[0]))
                    distances.append(int(parts[1]))

    return metadata, np.array(positions), np.array(distances)


def plot_single_gene_histogram(gene_name, histogram_file, gene_depth_info, output_file):
    """
    Create improved histogram plot for a single gene.
    Single panel with better binning and smoothed density.
    """
    metadata, positions, distances = read_histogram_data(histogram_file)

    # Create figure with single large panel
    fig, ax = plt.subplots(1, 1, figsize=(14, 7))

    # Adaptive binning based on gene length and read count
    gene_length = gene_depth_info['length_kb'] * 1000
    n_reads = len(distances)

    # Aim for ~20-30 bins regardless of gene length
    bin_size = max(50, int(gene_length / 25))  # At least 50bp bins
    bins = np.arange(distances.min() - bin_size, distances.max() + 2*bin_size, bin_size)

    # Plot histogram
    counts, edges, patches = ax.hist(distances, bins=bins, color='#3A86FF',
                                      alpha=0.6, edgecolor='black', linewidth=1.2,
                                      label=f'Observed (n={n_reads})')

    # Add smoothed density curve if we have enough data
    if SCIPY_AVAILABLE and n_reads >= 20:
        # Smooth the histogram
        bin_centers = (edges[:-1] + edges[1:]) / 2
        smoothed = gaussian_filter1d(counts, sigma=1.5)
        ax.plot(bin_centers, smoothed, color='#FB5607', linewidth=3,
                label='Smoothed trend', alpha=0.8)

    # Mark TSS
    ax.axvline(0, color='#D62828', linestyle='--', linewidth=2.5,
               label='TSS', alpha=0.8, zorder=10)

    # Add expected pattern annotation
    y_max = ax.get_ylim()[1]
    ax.axvspan(0, 100, alpha=0.15, color='green',
               label='Expected pausing region')

    # Styling
    ax.set_xlabel('Distance from TSS (bp)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Read Count per Bin', fontsize=14, fontweight='bold')
    ax.set_title(f"{gene_name} - PRO-seq 3' End Distribution\n"
                 f"Bin size: {bin_size} bp | Chr{metadata['chrom'].replace('chr', '')} | "
                 f"Strand: {metadata['strand']}",
                 fontsize=15, fontweight='bold', pad=20)
    ax.legend(fontsize=11, loc='upper right', framealpha=0.9)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Statistics box
    stats_text = (
        f"STATISTICS\n"
        f"==================\n"
        f"Total Reads: {n_reads}\n"
        f"Gene Length: {gene_depth_info['length_kb']:.2f} kb\n"
        f"Depth: {gene_depth_info['depth_per_kb']:.2f} reads/kb\n"
        f"Coordinates: {metadata['coords']}\n\n"
        f"WARNING: LOW COVERAGE\n"
        f"Noisy distribution due to\n"
        f"insufficient read depth.\n"
        f"Ideal: >1000 reads/gene"
    )

    ax.text(0.02, 0.97, stats_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow',
                     alpha=0.8, edgecolor='black', linewidth=1.5),
            family='monospace')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def plot_comparison_panel(gene_depth_df, output_file):
    """Multi-panel comparison with improved binning."""
    top_genes = gene_depth_df.head(TOP_N)

    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(4, 3, figure=fig, hspace=0.4, wspace=0.35)

    for idx, (_, row) in enumerate(top_genes.iterrows()):
        gene_name = row['gene_name']
        histogram_file = f"{HISTOGRAM_DIR}/{gene_name}_histogram.txt"

        if not os.path.exists(histogram_file):
            continue

        metadata, positions, distances = read_histogram_data(histogram_file)

        row_idx = idx // 3
        col_idx = idx % 3
        ax = fig.add_subplot(gs[row_idx, col_idx])

        # Adaptive binning
        gene_length = row['length_kb'] * 1000
        bin_size = max(100, int(gene_length / 15))
        bins = np.arange(distances.min() - bin_size, distances.max() + 2*bin_size, bin_size)

        # Plot
        ax.hist(distances, bins=bins, color='#3A86FF', alpha=0.7, edgecolor='black')
        ax.axvline(0, color='#D62828', linestyle='--', linewidth=2, alpha=0.8)
        ax.axvspan(0, 100, alpha=0.1, color='green')

        ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
        ax.set_ylabel('Read Count', fontsize=10)
        ax.set_title(f"#{idx+1}: {gene_name}\n{row['read_count']} reads | "
                     f"{row['depth_per_kb']:.1f} reads/kb",
                     fontsize=11, fontweight='bold')
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add warning if very low coverage
        if row['read_count'] < 100:
            ax.text(0.95, 0.95, 'Low\ncoverage', transform=ax.transAxes,
                    fontsize=8, ha='right', va='top',
                    bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    fig.suptitle('Top 10 Genes by PRO-seq Depth - 3\' End Distribution\n'
                 'WARNING: Chr22 data only | Low read depth causes noisy distributions',
                 fontsize=16, fontweight='bold', y=0.995)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def plot_depth_comparison(gene_depth_df, output_file):
    """Bar chart comparing depths."""
    top_genes = gene_depth_df.head(TOP_N)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    colors = plt.cm.viridis(np.linspace(0.3, 0.9, TOP_N))

    # Panel 1: Depth
    bars1 = ax1.barh(range(TOP_N), top_genes['depth_per_kb'], color=colors,
                      edgecolor='black', linewidth=1.5)
    ax1.set_yticks(range(TOP_N))
    ax1.set_yticklabels(top_genes['gene_name'], fontsize=11, fontweight='bold')
    ax1.set_xlabel('Sequencing Depth (reads/kb)', fontsize=12, fontweight='bold')
    ax1.set_title('Normalized Depth (accounts for gene length)',
                  fontsize=13, fontweight='bold', pad=15)
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    for i, (bar, val) in enumerate(zip(bars1, top_genes['depth_per_kb'])):
        ax1.text(val + 2, i, f'{val:.1f}', va='center', fontsize=10, fontweight='bold')

    # Panel 2: Total reads
    bars2 = ax2.barh(range(TOP_N), top_genes['read_count'], color=colors,
                      edgecolor='black', linewidth=1.5)
    ax2.set_yticks(range(TOP_N))
    ax2.set_yticklabels(top_genes['gene_name'], fontsize=11, fontweight='bold')
    ax2.set_xlabel('Total Read Count', fontsize=12, fontweight='bold')
    ax2.set_title('Absolute Read Count (not normalized)',
                  fontsize=13, fontweight='bold', pad=15)
    ax2.invert_yaxis()
    ax2.grid(axis='x', alpha=0.3, linestyle='--')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    for i, (bar, val) in enumerate(zip(bars2, top_genes['read_count'])):
        ax2.text(val + 5, i, f'{int(val)}', va='center', fontsize=10, fontweight='bold')

    plt.suptitle('Top 10 Genes - Sequencing Depth Comparison\n'
                 'WARNING: Low overall coverage (chr22 only, 0.94% mapping rate)',
                 fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def main():
    """Main visualization."""
    print("=" * 80)
    print("PRO-seq Visualization (Improved for Low-Coverage Data)")
    print("=" * 80)

    gene_depth_df = pd.read_csv(GENE_DEPTH_FILE, sep='\t')
    top_genes = gene_depth_df.head(TOP_N)

    print(f"\nTop {TOP_N} genes:")
    for idx, row in top_genes.iterrows():
        print(f"  {idx+1}. {row['gene_name']:15s} {row['depth_per_kb']:6.2f} reads/kb "
              f"({int(row['read_count'])} reads)")

    # Individual plots
    print(f"\n{'='*80}")
    print("Creating individual gene plots...")
    print('='*80)

    for idx, row in top_genes.iterrows():
        gene_name = row['gene_name']
        histogram_file = f"{HISTOGRAM_DIR}/{gene_name}_histogram.txt"

        if not os.path.exists(histogram_file):
            print(f"  WARNING: Skipping {gene_name} (file not found)")
            continue

        output_file = f"{OUTPUT_DIR}/{gene_name}_histogram_v2.png"
        gene_depth_info = {
            'length_kb': row['length_kb'],
            'depth_per_kb': row['depth_per_kb'],
            'read_count': row['read_count']
        }

        plot_single_gene_histogram(gene_name, histogram_file, gene_depth_info, output_file)

    # Comparison panel
    print(f"\n{'='*80}")
    print("Creating comparison panel...")
    print('='*80)
    comparison_file = f"{OUTPUT_DIR}/top10_comparison_v2.png"
    plot_comparison_panel(gene_depth_df, comparison_file)

    # Depth comparison
    print(f"\n{'='*80}")
    print("Creating depth comparison...")
    print('='*80)
    depth_file = f"{OUTPUT_DIR}/depth_comparison_v2.png"
    plot_depth_comparison(gene_depth_df, depth_file)

    print(f"\n{'='*80}")
    print("COMPLETE!")
    print('='*80)
    print(f"\nOutput: {OUTPUT_DIR}/")
    print(f"   - {TOP_N} individual gene histograms (*_v2.png)")
    print(f"   - 1 comparison panel (top10_comparison_v2.png)")
    print(f"   - 1 depth chart (depth_comparison_v2.png)")
    print("\nNOTE: Distributions are noisy due to low read coverage.")
    print("   Ideal PRO-seq requires >10M mapped reads (we have only 135k).")
    print('='*80)


if __name__ == "__main__":
    main()
