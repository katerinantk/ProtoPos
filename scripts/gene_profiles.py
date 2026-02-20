#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: gene_profiles.py.
Description:
    Plots single-base-resolution Pol II density profiles for individual genes
    from a 3'-end BAM file.  Each subplot shows the distribution of Pol II
    active-site positions along the gene body, with TSS and TES marked.

Dependencies:
    * Python >= 3.10, Python <= 3.13.
    * pysam, pandas, numpy, matplotlib

Usage:
    python scripts/gene_profiles.py
"""

import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from collections import defaultdict

# -------------------------------------------------------------------------
"""Inputs."""
# -------------------------------------------------------------------------

bam_path = "alignments/SRR11793827.3prime.sorted.bam"
gtf_path = "data/gencode.v47.annotation.gtf"
output_figure = "results/gene_profiles.png"

N_GENES = 10
SEED = 42


# -------------------------------------------------------------------------
"""Parse protein-coding genes from the GTF."""
# -------------------------------------------------------------------------

def parse_genes(gtf_path):
    """Parse one TSS per protein-coding gene (longest transcript).

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
            chrom = fields[0]
            if gid not in genes or length > genes[gid]["length"]:
                genes[gid] = {
                    "gene_name": name, "chrom": chrom, "strand": strand,
                    "start": start, "end": end, "length": length,
                }
    return pd.DataFrame(genes.values())


# -------------------------------------------------------------------------
"""Count Pol II reads per gene and select top candidates."""
# -------------------------------------------------------------------------

def count_reads_per_gene(bam_path, genes_df):
    """Count strand-concordant Pol II reads overlapping each gene.

    Args:
        bam_path (str): Path to indexed 3'-end BAM.
        genes_df (pd.DataFrame): Gene table with start, end, strand columns.

    Returns:
        pd.Series: Read counts indexed like genes_df.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    counts = []
    for _, row in genes_df.iterrows():
        n = 0
        for read in bam.fetch(row["chrom"], row["start"] - 1, row["end"]):
            # PRO-seq R1 is antisense to the RNA.
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
        # Strand filter: PRO-seq R1 is antisense to the RNA.
        if strand == "+" and not read.is_reverse:
            continue
        if strand == "-" and read.is_reverse:
            continue

        pos0 = read.reference_start
        idx = pos0 - (start - 1)
        if 0 <= idx < length:
            counts[idx] += 1

    bam.close()

    # Flip minus-strand genes so the plot reads 5' to 3'.
    if strand == "-":
        counts = counts[::-1]

    return counts


# -------------------------------------------------------------------------
"""Main."""
# -------------------------------------------------------------------------

def main():
    print("Parsing annotation...")
    genes = parse_genes(gtf_path)

    print("Counting reads per gene...")
    genes["read_count"] = count_reads_per_gene(bam_path, genes)

    # Select 10 random genes from the top 20 by read count.
    top20 = genes.sort_values("read_count", ascending=False).head(20)
    random.seed(SEED)
    selected = top20.sample(n=N_GENES, random_state=SEED).sort_values(
        "read_count", ascending=False
    )

    print(f"Selected {N_GENES} genes:")
    for _, row in selected.iterrows():
        print(f"  {row['gene_name']:15s}  ({row['strand']})  {row['read_count']:>4} reads  "
              f"{row['length']:>7,} bp")

    # Build the figure: 2 columns x 5 rows.
    fig, axes = plt.subplots(5, 2, figsize=(14, 16))
    axes = axes.flatten()

    smoothing_window = 201

    for i, (_, row) in enumerate(selected.iterrows()):
        ax = axes[i]
        profile = gene_profile(
            bam_path, row["chrom"], row["strand"],
            row["start"], row["end"],
        )

        # X-axis in kb from TSS.
        x_bp = np.arange(len(profile))
        x_kb = x_bp / 1000.0

        # Smooth for display.
        smooth = pd.Series(profile).rolling(
            smoothing_window, center=True, min_periods=1
        ).mean().values

        ax.fill_between(x_kb, smooth, alpha=0.3, color="#2c7bb6")
        ax.plot(x_kb, smooth, linewidth=0.8, color="#2c7bb6")

        # Mark TSS and TES.
        ax.axvline(0, color="#d7191c", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.axvline(x_kb[-1], color="#4dac26", linewidth=0.8, linestyle="--", alpha=0.6)

        ax.set_title(
            f"{row['gene_name']}  ({row['strand']} strand, "
            f"{row['read_count']} reads, {row['length']:,} bp)",
            fontsize=10, fontweight="bold",
        )
        ax.set_xlabel("Distance from TSS (kb)", fontsize=9)
        ax.set_ylabel("Pol II density", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=8)

    fig.suptitle(
        "Pol II density profiles  \u2014  qPRO-seq K562 (chr22)\n"
        "Red dashed = TSS, green dashed = TES",
        fontsize=13, fontweight="bold", y=0.995,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(output_figure, dpi=200)
    plt.close()
    print(f"\nWrote {output_figure}")


if __name__ == "__main__":
    main()
