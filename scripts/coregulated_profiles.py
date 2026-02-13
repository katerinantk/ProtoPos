#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: coregulated_profiles.py.
Description:
    Produces Pol II density figures for co-regulated genes in the 22q11.2
    region of chromosome 22.  Generates three panels:
      1. Individual gene profiles (subplots).
      2. Group metagene (all genes normalised and averaged).
      3. Top high-depth genes (>1000 reads) with dense Pol II signal.

Dependencies:
    * Python >= 3.10, Python <= 3.13.
    * pysam, pandas, numpy, matplotlib

Usage:
    python scripts/coregulated_profiles.py
"""

import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# -------------------------------------------------------------------------
"""Inputs."""
# -------------------------------------------------------------------------

bam_path = "alignments/SRR11793827.3prime.sorted.bam"
gtf_path = "data/gencode.v47.annotation.gtf"

# 22q11.2 region boundaries (approximate, GRCh38).
REGION_START = 18_900_000
REGION_END = 22_000_000

# Metagene window around TSS.
UPSTREAM = 500
DOWNSTREAM = 2000

# Output paths.
output_individual = "results/22q11_individual_profiles.png"
output_metagene = "results/22q11_group_metagene.png"
output_high_depth = "results/high_depth_genes.png"


# -------------------------------------------------------------------------
"""Parse protein-coding genes from the GTF."""
# -------------------------------------------------------------------------

def parse_genes(gtf_path):
    """Parse one record per protein-coding gene (longest transcript).

    Args:
        gtf_path (str): Path to GENCODE GTF.

    Returns:
        pd.DataFrame: Columns gene_name, chrom, strand, start, end, length.
    """
    genes = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = {}
            for tok in fields[8].strip().split(";"):
                tok = tok.strip()
                if not tok or " " not in tok:
                    continue
                k, v = tok.split(" ", 1)
                attrs[k] = v.strip('"')
            gt = attrs.get("transcript_type") or attrs.get("gene_type")
            if gt != "protein_coding":
                continue
            gid = attrs.get("gene_id")
            name = attrs.get("gene_name", gid)
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            length = end - start + 1
            if gid not in genes or length > genes[gid]["length"]:
                genes[gid] = {
                    "gene_name": name, "chrom": "chr22", "strand": strand,
                    "start": start, "end": end, "length": length,
                }
    return pd.DataFrame(genes.values())


# -------------------------------------------------------------------------
"""Count strand-concordant Pol II reads per gene."""
# -------------------------------------------------------------------------

def count_reads(bam_path, genes_df):
    """Count PRO-seq reads overlapping each gene (strand-aware).

    Args:
        bam_path (str): Path to indexed 3'-end BAM.
        genes_df (pd.DataFrame): Gene table.

    Returns:
        pd.Series: Read counts indexed like genes_df.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    counts = []
    for _, row in genes_df.iterrows():
        n = 0
        for read in bam.fetch("chr22", row["start"] - 1, row["end"]):
            if row["strand"] == "+" and not read.is_reverse:
                continue
            if row["strand"] == "-" and read.is_reverse:
                continue
            n += 1
        counts.append(n)
    bam.close()
    return pd.Series(counts, index=genes_df.index)


# -------------------------------------------------------------------------
"""Build a per-base profile for one gene."""
# -------------------------------------------------------------------------

def gene_profile(bam_path, chrom, strand, start, end):
    """Return per-base Pol II counts across the gene, oriented 5' to 3'.

    Args:
        bam_path (str): Path to indexed 3'-end BAM.
        chrom (str): Chromosome.
        strand (str): '+' or '-'.
        start (int): 1-based gene start (leftmost).
        end (int): 1-based gene end (rightmost).

    Returns:
        np.ndarray: Counts array, length = end - start + 1, oriented 5'->3'.
    """
    length = end - start + 1
    counts = np.zeros(length, dtype=float)
    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch(chrom, start - 1, end):
        if strand == "+" and not read.is_reverse:
            continue
        if strand == "-" and read.is_reverse:
            continue
        pos0 = read.reference_start
        idx = pos0 - (start - 1)
        if 0 <= idx < length:
            counts[idx] += 1

    bam.close()

    if strand == "-":
        counts = counts[::-1]

    return counts


# -------------------------------------------------------------------------
"""Build a metagene profile around the TSS."""
# -------------------------------------------------------------------------

def tss_profile(bam_path, chrom, strand, tss, upstream, downstream):
    """Count Pol II reads in a window around the TSS.

    Args:
        bam_path (str): Path to indexed 3'-end BAM.
        chrom (str): Chromosome.
        strand (str): '+' or '-'.
        tss (int): 1-based TSS coordinate.
        upstream (int): Bases upstream of TSS to include.
        downstream (int): Bases downstream of TSS to include.

    Returns:
        np.ndarray: Counts per position, length = upstream + downstream + 1.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    window_len = upstream + downstream + 1
    counts = np.zeros(window_len, dtype=float)
    tss0 = tss - 1

    if strand == "+":
        region_start = max(0, tss0 - upstream)
        region_end = tss0 + downstream + 1
    else:
        region_start = max(0, tss0 - downstream)
        region_end = tss0 + upstream + 1

    for read in bam.fetch(chrom, region_start, region_end):
        pos0 = read.reference_start
        if strand == "+":
            rel = pos0 - tss0
        else:
            rel = tss0 - pos0

        read_is_minus = read.is_reverse
        gene_is_minus = (strand == "-")
        if read_is_minus == gene_is_minus:
            continue

        idx = rel + upstream
        if 0 <= idx < window_len:
            counts[idx] += 1

    bam.close()
    return counts


# -------------------------------------------------------------------------
"""Figure 1: Individual profiles for 22q11.2 genes."""
# -------------------------------------------------------------------------

def plot_individual(genes_df, output_path):
    """Plot individual Pol II density profiles for 22q11.2 genes.

    Args:
        genes_df (pd.DataFrame): Gene table with read_count column.
        output_path (str): Path for the output PNG.
    """
    n = len(genes_df)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.5 * nrows))
    axes = axes.flatten()

    for i, (_, row) in enumerate(genes_df.iterrows()):
        ax = axes[i]
        profile = gene_profile(
            bam_path, row["chrom"], row["strand"],
            row["start"], row["end"],
        )
        x_kb = np.arange(len(profile)) / 1000.0

        smooth_win = max(51, min(301, len(profile) // 100 * 2 + 1))
        if smooth_win % 2 == 0:
            smooth_win += 1
        smooth = pd.Series(profile).rolling(
            smooth_win, center=True, min_periods=1
        ).mean().values

        ax.fill_between(x_kb, smooth, alpha=0.3, color="#2c7bb6")
        ax.plot(x_kb, smooth, linewidth=0.8, color="#2c7bb6")
        ax.axvline(0, color="#d7191c", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.axvline(x_kb[-1], color="#4dac26", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.set_title(
            f"{row['gene_name']}  ({row['strand']}, "
            f"{row['read_count']:,} reads, {row['length']:,} bp)",
            fontsize=9, fontweight="bold",
        )
        ax.set_xlabel("Distance from TSS (kb)", fontsize=8)
        ax.set_ylabel("Pol II density", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=7)

    # Hide unused axes.
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Pol II density \u2014 22q11.2 region genes  \u2014  qPRO-seq K562 (chr22)\n"
        "Red = TSS, green = TES",
        fontsize=13, fontweight="bold", y=1.0,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Wrote {output_path}")


# -------------------------------------------------------------------------
"""Figure 2: Group metagene for 22q11.2 genes."""
# -------------------------------------------------------------------------

def plot_group_metagene(genes_df, output_path):
    """Plot averaged metagene profile around the TSS for 22q11.2 genes.

    Args:
        genes_df (pd.DataFrame): Gene table.
        output_path (str): Path for the output PNG.
    """
    window_len = UPSTREAM + DOWNSTREAM + 1
    all_profiles = []

    for _, row in genes_df.iterrows():
        tss = row["start"] if row["strand"] == "+" else row["end"]
        profile = tss_profile(
            bam_path, row["chrom"], row["strand"],
            tss, UPSTREAM, DOWNSTREAM,
        )
        gene_body_mean = profile[UPSTREAM + 250: UPSTREAM + DOWNSTREAM].mean()
        if gene_body_mean > 0:
            all_profiles.append(profile / gene_body_mean)

    if not all_profiles:
        print("  No genes with signal for metagene.")
        return

    metagene = np.mean(all_profiles, axis=0)
    positions = np.arange(-UPSTREAM, DOWNSTREAM + 1)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    smooth_window = 31
    smooth = pd.Series(metagene).rolling(
        smooth_window, center=True, min_periods=1
    ).mean().values

    ax.fill_between(positions, smooth, alpha=0.25, color="#2c7bb6")
    ax.plot(positions, smooth, linewidth=1.2, color="#2c7bb6")
    ax.axvline(0, color="#d7191c", linewidth=1, linestyle="--", alpha=0.7, label="TSS")
    ax.axvspan(20, 60, alpha=0.08, color="#fdae61", label="Pause region (+20 to +60)")

    ax.set_xlabel("Distance from TSS (bp)", fontsize=12)
    ax.set_ylabel("Normalised Pol II density\n(fold over gene body)", fontsize=12)
    ax.set_title(
        f"22q11.2 region metagene ({len(all_profiles)} genes)"
        "  \u2014  qPRO-seq K562 (chr22)",
        fontsize=13, fontweight="bold",
    )
    ax.set_xlim(-UPSTREAM, DOWNSTREAM)
    ax.legend(fontsize=10, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda x, _: f"+{int(x)}" if x > 0 else str(int(x))
    ))

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"  Wrote {output_path}")


# -------------------------------------------------------------------------
"""Figure 3: High-depth genes (reads/kb > 20, >1000 total reads)."""
# -------------------------------------------------------------------------

def plot_high_depth(genes_df, output_path):
    """Plot profiles for genes with very high Pol II density.

    Args:
        genes_df (pd.DataFrame): Gene table with read_count and reads_per_kb.
        output_path (str): Path for the output PNG.
    """
    selection = genes_df[
        (genes_df["read_count"] >= 1000) & (genes_df["reads_per_kb"] >= 20)
    ].sort_values("reads_per_kb", ascending=False).head(9)

    if selection.empty:
        print("  No genes meet high-depth criteria.")
        return

    n = len(selection)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.5 * nrows))
    axes = axes.flatten()

    for i, (_, row) in enumerate(selection.iterrows()):
        ax = axes[i]
        profile = gene_profile(
            bam_path, row["chrom"], row["strand"],
            row["start"], row["end"],
        )
        x_kb = np.arange(len(profile)) / 1000.0

        smooth_win = max(51, min(201, len(profile) // 100 * 2 + 1))
        if smooth_win % 2 == 0:
            smooth_win += 1
        smooth = pd.Series(profile).rolling(
            smooth_win, center=True, min_periods=1
        ).mean().values

        ax.fill_between(x_kb, smooth, alpha=0.3, color="#d95f02")
        ax.plot(x_kb, smooth, linewidth=0.8, color="#d95f02")
        ax.axvline(0, color="#d7191c", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.axvline(x_kb[-1], color="#4dac26", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.set_title(
            f"{row['gene_name']}  ({row['strand']}, "
            f"{row['read_count']:,} reads, {row['reads_per_kb']:.0f} reads/kb)",
            fontsize=9, fontweight="bold",
        )
        ax.set_xlabel("Distance from TSS (kb)", fontsize=8)
        ax.set_ylabel("Pol II density", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=7)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "High-density Pol II profiles (\u22651000 reads, \u226520 reads/kb)"
        "  \u2014  qPRO-seq K562 (chr22)\n"
        "Red = TSS, green = TES",
        fontsize=13, fontweight="bold", y=1.0,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Wrote {output_path}")


# -------------------------------------------------------------------------
"""Main."""
# -------------------------------------------------------------------------

def main():
    print("Parsing annotation...")
    genes = parse_genes(gtf_path)
    print(f"  {len(genes)} protein-coding genes on chr22.")

    print("Counting reads per gene...")
    genes["read_count"] = count_reads(bam_path, genes)
    genes["reads_per_kb"] = genes["read_count"] / (genes["length"] / 1000.0)

    # 22q11.2 region genes with >1000 reads.
    region_genes = genes[
        (genes["start"] >= REGION_START) &
        (genes["end"] <= REGION_END) &
        (genes["read_count"] >= 1000)
    ].sort_values("read_count", ascending=False)

    print(f"\n22q11.2 genes with >1000 reads: {len(region_genes)}")
    for _, r in region_genes.iterrows():
        print(f"  {r['gene_name']:15s}  {r['strand']}  {r['read_count']:>5,} reads  "
              f"{r['length']:>8,} bp  {r['reads_per_kb']:.1f} reads/kb")

    print("\nFigure 1: Individual 22q11.2 gene profiles...")
    plot_individual(region_genes, output_individual)

    print("Figure 2: 22q11.2 group metagene...")
    plot_group_metagene(region_genes, output_metagene)

    print("Figure 3: High-depth genes...")
    plot_high_depth(genes, output_high_depth)

    print("\nDone.")


if __name__ == "__main__":
    main()
