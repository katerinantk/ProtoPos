#!/usr/bin/env python3

import pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Inputs
# -----------------------------
gtf_path = "data/gencode.v47.annotation.gtf"
bw_plus  = "data/GSE181161_PROseq_K562_hg38_pl.bigWig"
bw_minus = "data/GSE181161_PROseq_K562_hg38_mn.bigWig"

upperwindow = 200   # upstream of TSS
lowerwindow = 200   # downstream of TSS

# Outputs
output_file = "pol2signal.tsv"
output_plot = "pol2violin.png"


# ------------------------------------------------------------------------------
# Function: mRNAtranscripts
# One row per gene (longest transcript), protein_coding only
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

            # Parse attributes column into a dict
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

            # keep only mRNA / protein coding
            if gene_type not in ["protein_coding", "protein-coding", "mRNA"]:
                continue

            start = int(start)
            end   = int(end)
            length = end - start + 1

            # TSS depends on strand
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
                "length": length,
            })

    df = pd.DataFrame(rows)
    print(f"protein_coding-like transcripts: {len(df)}")

    # longest transcript per gene
    df_longest = (
        df.sort_values("length", ascending=False)
          .groupby("gene_id", as_index=False)
          .first()
    )

    print(f"longest transcript kept: {len(df_longest)}")
    return df_longest


# ------------------------------------------------------------------------------
# Function: tss_signal
# Get mean Pol2 signal in a window around the TSS
# ------------------------------------------------------------------------------
def tss_signal(bw_plus_handle, bw_minus_handle, chrom, strand, tss, up, down):
    # Choose window & bigWig depending on strand
    if strand == "+":
        start_bp = tss - up
        end_bp   = tss + down
        bw       = bw_plus_handle
    else:
        # for minus, "upstream" is higher coords
        start_bp = tss - down
        end_bp   = tss + up
        bw       = bw_minus_handle

    # convert to 0-based half-open
    start0 = max(0, start_bp - 1)
    end0   = max(0, end_bp)

    if start0 >= end0:
        return np.nan

    try:
        vals = bw.values(chrom, start0, end0, numpy=True)
    except RuntimeError:
        # chromosome not found in bigWig
        return np.nan

    if vals is None or len(vals) == 0:
        return np.nan

    vals = np.array(vals, dtype=float)

    # average Pol2 signal in this window
    return np.nanmean(vals)


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main():
    # annotation
    annotation = mRNAtranscripts(gtf_path)
    print("annotation done")

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
            down=lowerwindow,
        )
        signals.append(sig)

    # close bigWigs
    bw_pl.close()
    bw_mn.close()
    print("bigWigs closed")
    print("signals collected")

    # put signal into table
    annotation["pol2_signal"] = signals

    # drop genes with no signal
    annotation = annotation.dropna(subset=["pol2_signal"])

    # (optional) keep top 1000 by signal, sorted by gene name
    annotation = (
        annotation.sort_values("pol2_signal", ascending=False)
                  .head(1000)
                  .sort_values("gene_name")
                  .reset_index(drop=True)
    )

    # write table
    annotation.to_csv(output_file, sep="\t", index=False)
    print(f"wrote the table: {output_file}")

    # make the violin plot
    plt.figure(figsize=(6, 6))
    plt.title("RNAPII signal around the TSS")
    plt.ylabel("RNAPII signal")
    plt.violinplot(annotation["pol2_signal"].values, showmeans=True)
    plt.xticks([1], ["All genes"])
    plt.tight_layout()
    plt.savefig(output_plot, dpi=150)
    plt.close()
    print(f"the plot is done: {output_plot}")


if __name__ == "__main__":
    main()

