#!/usr/bin/env bash
# =============================================================================
# PRO-seq Gene Depth Analysis (Pure Bash)
# =============================================================================
# Authors: Katerina Ntouka and Christos Botos
# Affiliation: Institute of Molecular Biology and Biotechnology
# Contact: botoschristos@gmail.com
#
# Description:
#   Analyzes PRO-seq 3' end sequencing depth across protein-coding genes.
#   Uses only bash, awk, and standard Unix tools (no Python required).
#
# Usage:
#   bash scripts/analyze_gene_depth.sh
# =============================================================================

set -euo pipefail

# Configuration
GTF="data/gencode.v47.annotation.gtf"
SAM="alignments/SRR15304570.3prime.sam"
OUTPUT="results/gene_depth.tsv"
NUM_GENES=100

echo "========================================================================"
echo "PRO-seq Gene Depth Analysis"
echo "========================================================================"
echo "GTF file:    ${GTF}"
echo "SAM file:    ${SAM}"
echo "Output:      ${OUTPUT}"
echo "Num genes:   ${NUM_GENES}"
echo ""

# Create results directory
mkdir -p results

echo "Step 1: Extracting ${NUM_GENES} protein-coding genes from GTF..."

# Extract protein-coding genes from GTF and get coordinates
awk -F'\t' '
BEGIN { count = 0 }
$3 == "gene" && $9 ~ /gene_type "protein_coding"/ {
    # Extract gene_name from attributes
    match($9, /gene_name "([^"]+)"/, name)
    match($9, /gene_id "([^"]+)"/, id)

    # Calculate gene length in kb
    length_bp = $5 - $4 + 1
    length_kb = length_bp / 1000.0

    # Print: gene_name, chrom, start, end, strand, gene_id, length_kb
    printf "%s\t%s\t%d\t%d\t%s\t%s\t%.3f\n", name[1], $1, $4, $5, $7, id[1], length_kb

    count++
    if (count >= '"${NUM_GENES}"') exit
}' "${GTF}" > results/genes.tmp

ACTUAL_GENES=$(wc -l < results/genes.tmp)
echo "  Found ${ACTUAL_GENES} protein-coding genes"

echo ""
echo "Step 2: Counting reads in each gene from SAM file..."

# Header for output
echo -e "gene_name\tgene_id\tchrom\tstart\tend\tstrand\tlength_kb\tread_count\tdepth_per_kb" > "${OUTPUT}"

# For each gene, count overlapping reads in SAM file
while IFS=$'\t' read -r gene_name chrom start end strand gene_id length_kb; do
    # Count reads overlapping this gene
    # SAM format: column 3 is chromosome, column 4 is position
    read_count=$(awk -F'\t' -v chr="${chrom}" -v s="${start}" -v e="${end}" '
        $1 !~ /^@/ && $3 == chr && $4 >= s && $4 <= e { count++ }
        END { print count + 0 }
    ' "${SAM}")

    # Calculate depth (reads per kb)
    depth=$(awk -v r="${read_count}" -v l="${length_kb}" 'BEGIN { printf "%.2f", (l > 0 ? r / l : 0) }')

    # Output line
    echo -e "${gene_name}\t${gene_id}\t${chrom}\t${start}\t${end}\t${strand}\t${length_kb}\t${read_count}\t${depth}"

done < results/genes.tmp >> "${OUTPUT}"

echo "  Read counts completed"

echo ""
echo "Step 3: Sorting genes by sequencing depth..."

# Sort by depth (column 9) in descending order, keeping header
(head -n1 "${OUTPUT}" && tail -n +2 "${OUTPUT}" | sort -t$'\t' -k9 -rn) > results/gene_depth_sorted.tsv
mv results/gene_depth_sorted.tsv "${OUTPUT}"

echo "  Sorted by depth"

echo ""
echo "========================================================================"
echo "Top 10 Genes by Sequencing Depth"
echo "========================================================================"

# Display top 10 (skip header, take first 10 data lines)
echo -e "Rank\tGene\tReads\tLength(kb)\tDepth(reads/kb)"
tail -n +2 "${OUTPUT}" | head -10 | awk -F'\t' '{ printf "%d\t%s\t%d\t%.2f\t%.2f\n", NR, $1, $8, $7, $9 }'

echo ""
echo "========================================================================"
echo "Summary Statistics"
echo "========================================================================"

# Calculate stats
total_genes=$(tail -n +2 "${OUTPUT}" | wc -l)
genes_with_reads=$(tail -n +2 "${OUTPUT}" | awk -F'\t' '$8 > 0' | wc -l)
mean_depth=$(tail -n +2 "${OUTPUT}" | awk -F'\t' '{ sum += $9; count++ } END { printf "%.2f", sum/count }')
median_depth=$(tail -n +2 "${OUTPUT}" | awk -F'\t' '{ print $9 }' | sort -n | awk '{ a[NR] = $1 } END { if (NR % 2 == 1) print a[(NR+1)/2]; else print (a[NR/2] + a[NR/2+1])/2 }')
max_depth=$(tail -n +2 "${OUTPUT}" | awk -F'\t' '{ if ($9 > max) max = $9 } END { printf "%.2f", max }')

echo "Total genes analyzed:     ${total_genes}"
echo "Genes with reads:         ${genes_with_reads}"
echo "Mean depth (reads/kb):    ${mean_depth}"
echo "Median depth (reads/kb):  ${median_depth}"
echo "Max depth (reads/kb):     ${max_depth}"

echo ""
echo "Results saved to: ${OUTPUT}"
echo "========================================================================"

# Clean up
rm -f results/genes.tmp
