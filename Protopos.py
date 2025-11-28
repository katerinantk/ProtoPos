#!/usr/bin/env python3

import pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pysam

# -----------------------------
# Inputs
# -----------------------------

sam1_path = "star/SRR15304569.Aligned.sortedByCoord.out.sam"
sam2_path = "star/SRR15304570.Aligned.sortedByCoord.out.sam"
gtf_path = "data/gencode.v47.annotation.gtf"

def get3ends(in_path, out_path, gtf_path):
    # check if the sam files exist
    if not os.path.exists(sam1_path):
        raise FileNotFoundError(f"The first sam file was not found: {bam_path}")

    if not os.path.exists(sam2_path):
        raise FileNotFoundError(f"The second sam file was not found: {bam_path}")

    #check if gtf file exists
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")


# Open the input sam file 
infile = pysam.AlignmentFile(in_path, "r") 
# open the output sam file 
outfile = pysam.AlignmentFile(out_path, "w", header=infile.header)
















































    

   
