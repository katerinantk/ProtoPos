#!/usr/bin/env bash
# =============================================================================
# Generate Gene Histograms for Top Genes
# =============================================================================
# Authors: Katerina Ntouka and Christos Botos
# Affiliation: Institute of Molecular Biology and Biotechnology
# Contact: botoschristos@gmail.com
#
# Description:
#   Generates histogram data for the top N genes by sequencing depth.
#   For each gene, creates a file showing the distribution of 3' end positions
#   relative to the TSS (Transcription Start Site).
#
# Usage:
#   bash scripts/generate_gene_histograms.sh
# =============================================================================

set -euo pipefail

# Configuration
GTF="data/gencode.v47.annotation.gtf"
SAM="alignments/SRR15304570.3prime.sam"
GENE_DEPTH="results/gene_depth.tsv"
TOP_N=10
HISTOGRAM_DIR="results/histograms"

echo "========================================================================"
echo "PRO-seq Gene Histogram Generation"
echo "========================================================================"
echo "Gene depth file: ${GENE_DEPTH}"
echo "Top N genes:     ${TOP_N}"
echo "Output dir:      ${HISTOGRAM_DIR}"
echo ""

# Create output directory
mkdir -p "${HISTOGRAM_DIR}"

# Get top N gene names
echo "Extracting top ${TOP_N} genes..."
TOP_GENES=$(tail -n +2 "${GENE_DEPTH}" | head -${TOP_N} | cut -f1)

echo "Top genes:"
echo "${TOP_GENES}" | nl
echo ""

# Process each gene
for gene_name in ${TOP_GENES}; do
    echo "Processing gene: ${gene_name}"

    # Get gene coordinates from GTF
    gene_info=$(awk -F'\t' -v gene="${gene_name}" '
        $3 == "gene" && $9 ~ /gene_name "'"${gene_name}"'"/ {
            match($9, /gene_id "([^"]+)"/, id)
            printf "%s\t%d\t%d\t%s\t%s", $1, $4, $5, $7, id[1]
            exit
        }' "${GTF}")

    if [[ -z "${gene_info}" ]]; then
        echo "  WARNING: Gene ${gene_name} not found in GTF, skipping"
        continue
    fi

    # Parse gene coordinates
    chrom=$(echo "${gene_info}" | cut -f1)
    start=$(echo "${gene_info}" | cut -f2)
    end=$(echo "${gene_info}" | cut -f3)
    strand=$(echo "${gene_info}" | cut -f4)
    gene_id=$(echo "${gene_info}" | cut -f5)

    # Determine TSS based on strand
    if [[ "${strand}" == "+" ]]; then
        tss=${start}
    else
        tss=${end}
    fi

    echo "  Chrom: ${chrom}, TSS: ${tss}, Strand: ${strand}"

    # Extract reads overlapping this gene and calculate distance from TSS
    output_file="${HISTOGRAM_DIR}/${gene_name}_histogram.txt"

    # Header
    echo "# Gene: ${gene_name}" > "${output_file}"
    echo "# Gene ID: ${gene_id}" >> "${output_file}"
    echo "# Chromosome: ${chrom}" >> "${output_file}"
    echo "# TSS: ${tss}" >> "${output_file}"
    echo "# Strand: ${strand}" >> "${output_file}"
    echo "# Gene coordinates: ${start}-${end}" >> "${output_file}"
    echo "#" >> "${output_file}"
    echo "# position = genomic coordinate" >> "${output_file}"
    echo "# distance_from_TSS = position - TSS (positive = downstream, negative = upstream)" >> "${output_file}"
    echo "#" >> "${output_file}"
    echo -e "position\tdistance_from_TSS" >> "${output_file}"

    # Extract reads and calculate distance from TSS
    awk -F'\t' -v chr="${chrom}" -v s="${start}" -v e="${end}" -v tss="${tss}" -v str="${strand}" '
        $1 !~ /^@/ && $3 == chr && $4 >= s && $4 <= e {
            pos = $4
            if (str == "+") {
                dist = pos - tss
            } else {
                dist = tss - pos
            }
            printf "%d\t%d\n", pos, dist
        }
    ' "${SAM}" >> "${output_file}"

    read_count=$(tail -n +12 "${output_file}" | wc -l)
    echo "  Reads found: ${read_count}"

    # Create simple histogram summary (bins of 100 bp)
    summary_file="${HISTOGRAM_DIR}/${gene_name}_summary.txt"
    echo "# Histogram summary for ${gene_name}" > "${summary_file}"
    echo "# Bin size: 100 bp" >> "${summary_file}"
    echo -e "bin_start\tbin_end\tcount" >> "${summary_file}"

    tail -n +12 "${output_file}" | awk -F'\t' '
        {
            dist = $2
            bin = int(dist / 100) * 100
            bins[bin]++
        }
        END {
            for (bin in bins) {
                printf "%d\t%d\t%d\n", bin, bin + 100, bins[bin]
            }
        }
    ' | sort -n -k1 >> "${summary_file}"

    echo "  Saved to: ${output_file}"
    echo "  Summary:  ${summary_file}"
    echo ""
done

echo "========================================================================"
echo "Histogram generation complete!"
echo "Output directory: ${HISTOGRAM_DIR}"
echo "========================================================================"
