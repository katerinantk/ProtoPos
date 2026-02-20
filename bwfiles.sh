#!/usr/bin/env bash
set -euo pipefail

BAM="${1:-alignments/SRR28248970.sorted.3prime.sorted.bam}"
CHROMS="${2:-data/genome.chrom.sizes}"
OUT="${3:-SRR28248970}"

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM"
    exit 1
fi

if [[ ! -f "$CHROMS" ]]; then
    echo "ERROR: Chrom sizes file not found: $CHROMS"
    exit 1
fi

# plus strand
bedtools genomecov -ibam "$BAM" -bg -strand + -g "$CHROMS" \
  | sort -k1,1 -k2,2n > "${OUT}.plus.bedGraph"

# minus strand
bedtools genomecov -ibam "$BAM" -bg -strand - -g "$CHROMS" \
  | sort -k1,1 -k2,2n > "${OUT}.minus.bedGraph"


bedGraphToBigWig "${OUT}.plus.bedGraph"  "$CHROMS" "${OUT}.plus.bw"
bedGraphToBigWig "${OUT}.minus.bedGraph" "$CHROMS" "${OUT}.minus.bw"


