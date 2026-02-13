#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: steady_state_figure.py.
Description:
    Produces a metagene plot of RNA Polymerase II density around the TSS from
    a 3'-end BAM file.  The figure shows the characteristic promoter-proximal
    pausing peak at +20 to +60 bp downstream of the TSS, which is a hallmark
    of steady-state PRO-seq signal.

Dependencies:
    * Python >= 3.10, Python <= 3.13.
    * pysam, pandas, numpy, matplotlib

Usage:
    python scripts/steady_state_figure.py
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

# Window around TSS for the metagene profile.
UPSTREAM = 500
DOWNSTREAM = 2000

# Minimum number of reads in the gene body to include a gene.
MIN_GENE_BODY_READS = 1

# Output paths.
output_figure = "results/steady_state_metagene.png"
output_table = "results/steady_state_metagene.tsv"


# -------------------------------------------------------------------------
"""Parse protein-coding TSSs from the GTF (one per gene, longest transcript)."""
# -------------------------------------------------------------------------

def parse_tss(gtf_path):
    """Parse one TSS per protein-coding gene from the GTF.

    Keeps the longest transcript per gene to obtain a single representative
    TSS.

    Args:
        gtf_path (str): Path to the GENCODE GTF file.

    Returns:
        pd.DataFrame: Columns chrom, strand, tss, gene_name.
    """
    rows = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            attrs = {}
            for token in fields[8].strip().split(";"):
                token = token.strip()
                if not token or " " not in token:
                    continue
                key, val = token.split(" ", 1)
                attrs[key] = val.strip('"')

            gene_type = (
                attrs.get("transcript_type")
                or attrs.get("gene_type")
            )
            if gene_type != "protein_coding":
                continue

            chrom = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            tss = start if strand == "+" else end
            length = end - start + 1

            rows.append({
                "gene_id": attrs.get("gene_id"),
                "gene_name": attrs.get("gene_name", ""),
                "chrom": chrom,
                "strand": strand,
                "tss": tss,
                "length": length,
            })

    df = pd.DataFrame(rows)
    # Keep the longest transcript per gene.
    df = (
        df.sort_values("length", ascending=False)
          .groupby("gene_id", as_index=False)
          .first()
    )
    return df


# -------------------------------------------------------------------------
"""Build per-base Pol II density vector from the 3'-end BAM."""
# -------------------------------------------------------------------------

def build_profile(bam_path, chrom, strand, tss, upstream, downstream):
    """Count Pol II reads in a window around the TSS.

    The window is oriented so that negative positions are upstream of the TSS
    and positive positions are downstream, regardless of strand.

    Args:
        bam_path (str): Path to the indexed 3'-end BAM.
        chrom (str): Chromosome name.
        strand (str): '+' or '-'.
        tss (int): 1-based TSS coordinate.
        upstream (int): Bases upstream of the TSS.
        downstream (int): Bases downstream of the TSS.

    Returns:
        np.ndarray: Read counts per position, length = upstream + downstream + 1.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    window_len = upstream + downstream + 1
    counts = np.zeros(window_len, dtype=float)

    # Convert 1-based TSS to 0-based for pysam.
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

        # Strand concordance for PRO-seq: R1 maps antisense to the RNA, so
        # a + strand gene produces reverse-strand reads and vice versa.
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
"""Main."""
# -------------------------------------------------------------------------

def main():
    print("Parsing annotation...")
    genes = parse_tss(gtf_path)
    print(f"  {len(genes)} protein-coding genes.")

    print("Building metagene profile...")
    window_len = UPSTREAM + DOWNSTREAM + 1
    all_profiles = []

    for _, row in genes.iterrows():
        profile = build_profile(
            bam_path, row["chrom"], row["strand"],
            row["tss"], UPSTREAM, DOWNSTREAM,
        )
        # Only include genes with enough signal in the gene body (+250 to +2000).
        gene_body_start = UPSTREAM + 250
        gene_body_end = UPSTREAM + DOWNSTREAM
        gene_body_signal = profile[gene_body_start:gene_body_end].sum()
        if gene_body_signal >= MIN_GENE_BODY_READS:
            all_profiles.append(profile)

    print(f"  {len(all_profiles)} genes with sufficient signal.")

    if not all_profiles:
        print("ERROR: No genes with sufficient signal. Check BAM or annotation.")
        return

    # Normalise each gene to its own gene body mean, then average across genes.
    normed = []
    for p in all_profiles:
        gene_body_mean = p[UPSTREAM + 250 : UPSTREAM + DOWNSTREAM].mean()
        if gene_body_mean > 0:
            normed.append(p / gene_body_mean)

    metagene = np.mean(normed, axis=0)
    positions = np.arange(-UPSTREAM, DOWNSTREAM + 1)

    # Write table.
    table = pd.DataFrame({"position": positions, "signal": metagene})
    table.to_csv(output_table, sep="\t", index=False)
    print(f"  Wrote {output_table}")

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(10, 4.5))

    # Smooth with a rolling window for a clean line.
    smooth_window = 51 if len(all_profiles) < 50 else 15
    smooth = pd.Series(metagene).rolling(smooth_window, center=True, min_periods=1).mean().values

    ax.fill_between(positions, smooth, alpha=0.25, color="#2c7bb6")
    ax.plot(positions, smooth, linewidth=1.2, color="#2c7bb6")

    # Mark TSS.
    ax.axvline(0, color="#d7191c", linewidth=1, linestyle="--", alpha=0.7, label="TSS")

    # Mark typical pause region.
    ax.axvspan(20, 60, alpha=0.08, color="#fdae61", label="Pause region (+20 to +60)")

    ax.set_xlabel("Distance from TSS (bp)", fontsize=12)
    ax.set_ylabel("Normalised Pol II density\n(fold over gene body)", fontsize=12)
    ax.set_title(
        "Steady-state RNAPII metagene profile  —  qPRO-seq K562 (chr22)",
        fontsize=13, fontweight="bold",
    )
    ax.set_xlim(-UPSTREAM, DOWNSTREAM)
    ax.legend(fontsize=10, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Format x-axis with +/- labels.
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda x, _: f"+{int(x)}" if x > 0 else str(int(x))
    ))

    plt.tight_layout()
    plt.savefig(output_figure, dpi=200)
    plt.close()
    print(f"  Wrote {output_figure}")
    print("Done.")


if __name__ == "__main__":
    main()
