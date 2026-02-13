#!/usr/bin/env python3
"""
Authors: Katerina Ntouka and Christos Botos.
Affiliation: Institute of Molecular Biology and Biotechnology.
Contact: botoschristos@gmail.com

Script Name: Protopos.py.
Description:
    Extracts the RNA Polymerase II active-site position from aligned PRO-seq
    reads.  For each mapped read the script creates a single-base (1M)
    alignment at the 5' end of the sequencing read, which corresponds to the
    3' end of the nascent RNA and therefore the Pol II active site.

Dependencies:
    • Python >= 3.10, Python <= 3.13.
    • pysam

Usage:
    python Protopos.py
"""

import pysam

# ----------------------------------------------------
"""Input BAM files produced by STAR."""
# ----------------------------------------------------

bam_files = [
    "star/SRR15304569.Aligned.sortedByCoord.out.bam",
    "star/SRR15304570.Aligned.sortedByCoord.out.bam",
]


def pol2_position(read):
    """Return the 0-based genomic position of the Pol II active site.

    In standard PRO-seq the sequencing read starts from the 3' end of the
    nascent RNA (adjacent to the 3' adapter).  The 5' end of the aligned
    read therefore marks the RNA 3' terminus, i.e. the Pol II active site.

    Args:
        read (pysam.AlignedSegment): A mapped read.

    Returns:
        int: 0-based coordinate of the polymerase position.
    """
    if read.is_reverse:
        # Minus-strand alignment: 5' end of the read is at the rightmost position.
        return read.reference_end - 1
    else:
        # Plus-strand alignment: 5' end of the read is at the leftmost position.
        return read.reference_start


# ----------------------------------------------------
"""Process each BAM file."""
# ----------------------------------------------------

for input_bam in bam_files:

    output_bam = input_bam.replace(
        ".Aligned.sortedByCoord.out.bam", ".3prime.bam"
    )

    print(f"Processing {input_bam} ...")

    infile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", header=infile.header)

    for read in infile:
        if read.is_unmapped:
            continue

        pos = pol2_position(read)

        new = pysam.AlignedSegment(outfile.header)
        new.query_name = read.query_name
        new.flag = read.flag
        new.reference_id = read.reference_id
        new.reference_start = pos
        new.cigarstring = "1M"
        new.mapping_quality = read.mapping_quality

        new.query_sequence = "N"
        new.query_qualities = pysam.qualitystring_to_array("!")

        outfile.write(new)

    infile.close()
    outfile.close()

    print(f"Finished {output_bam}")
