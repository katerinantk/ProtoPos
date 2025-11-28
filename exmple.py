usr/bin/env python3

import pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pysam

# ------------------------------
# input file 
# ------------------------------

INPUT_SAM = "input.sam"
OUTPUT_SAM = "output_3prime.sam"
# ------------------------------

infile = pysam.AlignmentFile(INPUT_SAM, "r")
outfile = pysam.AlignmentFile(OUTPUT_SAM, "wh", header=infile.header)

def three_prime_position(read):
    
    #Compute 3' end of the read.
    
    if read.is_reverse:
        return read.reference_start
    else:
        return read.reference_end - 1

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

    # minimal placeholder sequence + quality
    new.query_sequence = "N"
    new.query_qualities = pysam.qualitystring_to_array("!")

    outfile.write(new)

infile.close()
outfile.close()

print("Done. Wrote:", OUTPUT_SAM)
