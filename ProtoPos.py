#!/usr/bin/env python3

import pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# These are the inputs
gtf_path = "data/gencode.v47.annotation.gtf"
bw_plus  = "data/GSE181161_PROseq_K562_hg38_pl.bigWig"
bw_minus = "data/GSE181161_PROseq_K562_hg38_mn.bigWig"

upperwindow = 200   # upstream of TSS
lowerwindow = 200   # downstream of TSS

# These are the outputs
output_file = "pol2signal.tsv"
output_plot = "pol2violin.png"


# ------------------------------------------------------------------------------
# Function: mRNAtranscripts
# This returns a dataframe with one row per gene (longest transcript)
# columns: gene_id, gene_name, chrom, strand, tss, length, etc.
# Only protein_coding / mRNA-like transcripts are kept.
# ------------------------------------------------------------------------------

def mRNAtranscripts(gtf_path):
    rows = []

    with open(gtf_path) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields

            if feature != "transcript":
                continue

            # This takes the attributes column from the gtf file and creates a dictionary
            attributes = {}
            for k in attrs.strip().split(";"):
                k = k.strip()
                if not k:
                    continue
                if " " in k:
                    key, val = k.split(" ", 1)
                    val = val.strip().strip('"')
                    attributes[key] = val

            gene_id       = attributes.get("gene_id")
            gene_name     = attributes.get("gene_name", gene_id)
            transcript_id = attributes.get("transcript_id")

            gene_type = (
                attributes.get("transcript_type")
                or attributes.get("transcript_biotype")
                or attributes.get("gene_type")
                or attributes.get("gene_biotype")
            )

            # keep only mRNA transcripts
            if gene_type not in ["protein_coding", "protein-coding", "mRNA"]:
                continue

            start = int(start)
            end   = int(end)
            length = end - start + 1  

            # find the TSS depending on the strand
            if strand == "+":
                tss = start
            else:
                tss = end

            rows.append({
                "gene_id": gene_id,
                "gene_name": gene_name,
                "transcript_id": transcript_id,
                "chrom": chrom,
                "strand": strand,
                "tss": tss,
                "tx_start": start,
                "tx_end": end,
                "length": length
            })

    df = pd.DataFrame(rows)

    print(f"protein_coding-like transcripts: {len(df)}")

    # this picks the longest transcript per gene
    df_longest = (
        df.sort_values("length", ascending=False)
          .groupby("gene_id", as_index=False)
          .first()
    )

    print(f"longest transcript kept: {len(df_longest)} ")

    return df_longest


# ------------------------------------------------------------------------------
# Function: tss_signal
# Get mean PRO-seq signal around TSS for a single gene.
# For + strand genes: use plus bigWig.
# For - strand genes: use minus bigWig.
# ------------------------------------------------------------------------------

def tss_signal(bw_plus_handle, bw_minus_handle, chrom, strand, tss, up, down):

    # define window around TSS 
    if strand == "+":
        start_bp = tss - up
        end_bp   = tss + down
        bw       = bw_plus_handle
    else:
        # minus strand, flip upstream,downstream
        start_bp = tss - down
        end_bp   = tss + up
        bw       = bw_minus_handle

    # convert 1-based inclusive coordinates to 0-based half-open for the bigWig files
    start0 = max(0, start_bp - 1)
    end0   = max(0, end_bp)

    try:
        vals = bw.values(chrom, start0, end0, numpy=True)
    except RuntimeError:
        # in case chromosome isnot found 
        return np.nan

    if vals is None or len(vals) == 0:
        return np.nan

    # average Pol2 signal in this window
    return np.nanmean(vals)


# ------------------------------------------------------------------------------
# Function: main
# ------------------------------------------------------------------------------

def main():
    # annotation 
    annotation = mRNAtranscripts(gtf_path)

    # open bigWigs
    bw_pl = pyBigWig.open(bw_plus)
    bw_mn = pyBigWig.open(bw_minus)

    # collect signal per gene
    signals = []
    for idx, row in annotation.iterrows():
        sig = tss_signal(
            bw_plus_handle=bw_pl,
            bw_minus_handle=bw_mn,
            chrom=row["chrom"],
            strand=row["strand"],
            tss=row["tss"],
            up=upperwindow,
            down=lowerwindow
        )
        signals.append(sig)

        

    # close bigWigs
    bw_pl.close()
    bw_mn.close()

    # put signal into table
    annotation["pol2_signal"] = signals

    # write table
    annotation.to_csv(output_file, sep="\t", index=False)
    print(f"[info] wrote table: {output_file}")

    # make violin plot
    plt.figure(figsize=(4,6))
    plt.violinplot(
        dataset=[annotation["pol2_signal"].dropna()],
        showmeans=True,
        showextrema=True,
        showmedians=True
    )
    plt.ylabel(f"mean PRO-seq signal at TSS (±{upperwindow} bp)")
    plt.title("K562 Pol2 promoter signal")
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    

    print("done")


if __name__ == "__main__":
    main()

