#!/usr/bin/env python3

import pysam

# ----------------------------------------------------
# input sam files
# ----------------------------------------------------

sam_files = [
    "star/SRR15304569.Aligned.sortedByCoord.out.sam",
    "star/SRR15304570.Aligned.sortedByCoord.out.sam"
]
  #compute 3' end of the read
def three_prime_position(read):
  
    if read.is_reverse:
        return read.reference_start
    else:
        return read.reference_end - 1

# process each sam file
for input_sam in sam_files:

    output_sam = input_sam.replace(".Aligned.sortedByCoord.out.sam", ".3prime.sam")

    print("processing the sam files")

    infile = pysam.AlignmentFile(input_sam, "r")
    outfile = pysam.AlignmentFile(output_sam, "wh", header=infile.header)

    for read in infile:
        if read.is_unmapped:
            continue

        tp = three_prime_position(read)

        new = pysam.AlignedSegment(outfile.header)
        new.query_name = read.query_name
        new.flag = read.flag
        new.reference_id = read.reference_id
        new.reference_start = tp
        new.cigarstring = "1M"
        new.mapping_quality = read.mapping_quality

        new.query_sequence = "n"
        new.query_qualities = pysam.qualitystring_to_array("!")
        

        outfile.write(new)

    infile.close()
    outfile.close()

    print(f"finished {output_sam}")


