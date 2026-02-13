#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: Protopos.py.
Description:
    Extracts the RNA Polymerase II active-site position from aligned
    nascent-RNA sequencing reads.  Supports both PRO-seq and GRO-seq
    library conventions:

    PRO-seq (antisense library):
        R1 maps antisense to the nascent RNA.  The 5' end of the aligned
        read corresponds to the 3' end of the RNA, i.e. the Pol II active
        site.
        - Forward alignment: active site = reference_start.
        - Reverse alignment: active site = reference_end - 1.

    GRO-seq (sense library):
        Reads map in the same direction as transcription.  The 3' end of
        the aligned read marks the polymerase active site.
        - Forward alignment: active site = reference_end - 1.
        - Reverse alignment: active site = reference_start.

    For each mapped read the script creates a single-base (1M) alignment
    at the active-site position.

Dependencies:
    * Python >= 3.10, Python <= 3.13.
    * pysam

Usage:
    python Protopos.py
"""

import os

import pysam


# -------------------------------------------------------------------------
"""Pol II active-site extraction."""
# -------------------------------------------------------------------------

def pol2_position(read, protocol="pro_seq"):
    """Return the 0-based genomic position of the Pol II active site.

    Args:
        read (pysam.AlignedSegment): A mapped read.
        protocol (str): Library protocol.
            ``"pro_seq"`` — PRO-seq / qPRO-seq (antisense R1).  The 5'
            end of the sequencing read marks the active site.
            ``"gro_seq"`` — GRO-seq (sense reads).  The 3' end of the
            sequencing read marks the active site.

    Returns:
        int: 0-based coordinate of the polymerase position.

    Raises:
        ValueError: If *protocol* is not recognised.
    """
    if protocol == "pro_seq":
        # PRO-seq: 5' end of the read = 3' of nascent RNA = active site.
        if read.is_reverse:
            return read.reference_end - 1
        else:
            return read.reference_start
    elif protocol == "gro_seq":
        # GRO-seq: 3' end of the read = active site (sense orientation).
        if read.is_reverse:
            return read.reference_start
        else:
            return read.reference_end - 1
    else:
        raise ValueError(
            f"Unknown protocol '{protocol}'. Use 'pro_seq' or 'gro_seq'."
        )


def extract_active_sites(input_bam, output_bam, protocol="pro_seq"):
    """Extract Pol II active-site positions into a 1 bp BAM file.

    Creates a new BAM where every mapped read is replaced by a single-base
    alignment at the polymerase active-site position determined by
    *protocol*.  The output BAM is sorted and indexed in place.

    Args:
        input_bam (str): Path to input sorted BAM file.
        output_bam (str): Path for the output 1 bp BAM file.  After
            sorting the file is written to this path (sorted and indexed).
        protocol (str): ``"pro_seq"`` or ``"gro_seq"``.

    Returns:
        dict: Statistics with keys ``'total_reads'``, ``'plus_strand'``,
            ``'minus_strand'``.
    """
    stats = {"total_reads": 0, "plus_strand": 0, "minus_strand": 0}

    infile = pysam.AlignmentFile(input_bam, "rb")
    tmp_bam = output_bam + ".tmp.bam"
    outfile = pysam.AlignmentFile(tmp_bam, "wb", header=infile.header)

    for read in infile:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        stats["total_reads"] += 1
        if read.is_reverse:
            stats["minus_strand"] += 1
        else:
            stats["plus_strand"] += 1

        pos = pol2_position(read, protocol)

        new = pysam.AlignedSegment(outfile.header)
        new.query_name = read.query_name
        new.flag = read.flag
        new.reference_id = read.reference_id
        new.reference_start = pos
        new.cigarstring = "1M"
        new.mapping_quality = read.mapping_quality
        new.query_sequence = "N"
        new.query_qualities = [0]

        outfile.write(new)

    infile.close()
    outfile.close()

    # Sort and index.
    pysam.sort("-o", output_bam, tmp_bam)
    pysam.index(output_bam)
    os.remove(tmp_bam)

    return stats


# -------------------------------------------------------------------------
"""Default BAM file list (edit to match your experiment)."""
# -------------------------------------------------------------------------

BAM_FILES = [
    "alignments/SRR15304569.dedup.bam",
    "alignments/SRR15304570.dedup.bam",
]


# -------------------------------------------------------------------------
"""Main."""
# -------------------------------------------------------------------------

def main():
    """Run active-site extraction on all BAM files in BAM_FILES."""
    for input_bam in BAM_FILES:
        output_bam = input_bam.replace(".dedup.bam", ".3prime.sorted.bam")
        print(f"Processing {input_bam} ...")
        stats = extract_active_sites(input_bam, output_bam, protocol="pro_seq")
        print(f"  {stats['total_reads']} reads "
              f"({stats['plus_strand']} +, {stats['minus_strand']} -)")
        print(f"  Wrote {output_bam}")


if __name__ == "__main__":
    main()
