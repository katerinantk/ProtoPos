#!/usr/bin/env bash
BAM="alignments/SRR28248970.sorted.3prime.sorted.bam"
CHROMS="data/genome.chrom.sizes"   
OUT="SRR28248970"

# plus strand
bedtools genomecov -ibam "$BAM" -bg -strand + -g "$CHROMS" \
  | sort -k1,1 -k2,2n > "${OUT}.plus.bedGraph"

# minus strand
bedtools genomecov -ibam "$BAM" -bg -strand - -g "$CHROMS" \
  | sort -k1,1 -k2,2n > "${OUT}.minus.bedGraph"


bedGraphToBigWig "${OUT}.plus.bedGraph"  "$CHROMS" "${OUT}.plus.bw"
bedGraphToBigWig "${OUT}.minus.bedGraph" "$CHROMS" "${OUT}.minus.bw"


