#!/usr/bin/env bash

# -----------------------------
# Inputs
# -----------------------------
SRA_IDS=("SRR15304569" "SRR15304570")
THREADS=8

DATA_DIR="data"
FASTQ_DIR="fastq"
FASTP_DIR="fastp"
QC_DIR="qc"
STAR_DIR="star"
STAR_INDEX_DIR="starindex"
BED_DIR="bed"
ENDS_DIR="ends"
BIGWIG_DIR="bigwig"

GENOME_FA="${DATA_DIR}/genome.fa"
GTF="${DATA_DIR}/gencode.v47.annotation.gtf"

# -----------------------------
# Check if the FASTA and GTF files exist
# -----------------------------
if [[ ! -f "${GENOME_FA}" ]]; then
    echo "ERROR: Genome file not found: ${GENOME_FA}"
    exit 1
fi

if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF file not found: ${GTF}"
    exit 1
fi

# -----------------------------
# Create output directories 
# -----------------------------
mkdir -p \
  "${FASTQ_DIR}" \
  "${FASTP_DIR}" \
  "${QC_DIR}" \
  "${STAR_DIR}" \
  "${STAR_INDEX_DIR}" \
  "${BED_DIR}" \
  "${ENDS_DIR}" \
  "${BIGWIG_DIR}"

# -----------------------------
# Build STAR genome index 
# -----------------------------
if [[ ! -f "${STAR_INDEX_DIR}/Genome" ]]; then
    STAR \
      --runThreadN "${THREADS}" \
      --runMode genomeGenerate \
      --genomeDir "${STAR_INDEX_DIR}" \
      --genomeFastaFiles "${GENOME_FA}" \
      --sjdbGTFfile "${GTF}" \
      --sjdbOverhang 75
fi

# -----------------------------
# Download / Fastqc / fastp / star 
# -----------------------------
for SRA in "${SRA_IDS[@]}"; do

    RAW="${FASTQ_DIR}/${SRA}.fastq.gz"
    TRIM="${FASTP_DIR}/${SRA}.trimmed.fastq.gz"

    HTML_REPORT="${QC_DIR}/${SRA}.fastp.html"
    JSON_REPORT="${QC_DIR}/${SRA}.fastp.json"

    # -------------------------
    # Download fastq
    # -------------------------
    if [[ ! -f "${RAW}" ]]; then
        fasterq-dump "${SRA}" -O "${FASTQ_DIR}" --threads "${THREADS}"
        gzip -f "${FASTQ_DIR}/${SRA}.fastq"
    else
        echo "Skipping download for ${SRA}: ${RAW} already exists."
    fi

    # -------------------------
    # Fastqc on raw reads
    # -------------------------
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${RAW}"

    # -------------------------
    # Trimming with fastp 
    # -------------------------
    fastp \
      -i "${RAW}" \
      -o "${TRIM}" \
      -w "${THREADS}" \
      --html "${HTML_REPORT}" \
      --json "${JSON_REPORT}"

    # -------------------------
    # Fastqc on trimmed reads
    # -------------------------
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${TRIM}"

    # -------------------------
    # star alignment 
    # -------------------------
    PREFIX="${STAR_DIR}/${SRA}."
    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX_DIR}" \
      --readFilesIn "${TRIM}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${PREFIX}" \
      --outSAMtype BAM SortedByCoordinate

    
done

# -----------------------------
# rRNA and tRNA BED regions 
# -----------------------------
RRNA_TRNA_BED="${DATA_DIR}/rRNAtRNA.bed"
RRNA_TRNA_UNIQ="${DATA_DIR}/rRNAtRNA.unique.bed"
RRNA_TRNA_MERGED="${DATA_DIR}/rRNAtRNA.merged.bed"

if [[ ! -f "${RRNA_TRNA_MERGED}" ]]; then

    awk '
        $3 == "exon" && $0 ~ /gene_type "rRNA"/ {print $1"\t"$4-1"\t"$5}
        $3 == "exon" && $0 ~ /gene_type "tRNA"/ {print $1"\t"$4-1"\t"$5}
    ' "${GTF}" > "${RRNA_TRNA_BED}"

    sort -u "${RRNA_TRNA_BED}" > "${RRNA_TRNA_UNIQ}"

    sort -k1,1 -k2,2n "${RRNA_TRNA_UNIQ}" \
      | bedtools merge -i - > "${RRNA_TRNA_MERGED}"
fi

# -----------------------------
# Chromosome sizes 
# -----------------------------
CHROM_SIZES="${DATA_DIR}/genome.chrom.sizes"

if [[ ! -f "${CHROM_SIZES}" ]]; then
    samtools faidx "${GENOME_FA}"
    cut -f1,2 "${GENOME_FA}.fai" > "${CHROM_SIZES}"
fi

# -----------------------------
# From mapped bams to bed to 3' ends 
# we filter rRNA and tRNA 
# -----------------------------
for SRA in "${SRA_IDS[@]}"; do

    # Input from bams from STAR
    BAM_SORTED="${STAR_DIR}/${SRA}.Aligned.sortedByCoord.out.bam"

    # Keep only mapped reads
    BAM_MAPPED="${STAR_DIR}/${SRA}.mapped.bam"
    samtools view -b -F 4 "${BAM_SORTED}" > "${BAM_MAPPED}"
    samtools index "${BAM_MAPPED}"

    # bam to bed 
    BED_ALL="${BED_DIR}/${SRA}.bed"
    bedtools bamtobed -i "${BAM_MAPPED}" > "${BED_ALL}"

    # get the 3' ends in a strand specific way 
    ENDS_3="${ENDS_DIR}/${SRA}.3end.bed"
    awk '
        ($6 == "+") {print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}
        ($6 == "-") {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6}
    ' "${BED_ALL}" > "${ENDS_3}"

    # Remove rRNA and tRNA regions
    ENDS_FILTERED="${ENDS_DIR}/${SRA}.filtered.3end.bed"
    bedtools intersect -v \
      -a "${ENDS_3}" \
      -b "${RRNA_TRNA_MERGED}" \
      > "${ENDS_FILTERED}"

    # Split by strand
    ENDS_PLUS="${ENDS_DIR}/${SRA}.filtered.plus.3end.bed"
    ENDS_MINUS="${ENDS_DIR}/${SRA}.filtered.minus.3end.bed"

    awk '$6 == "+"' "${ENDS_FILTERED}"  > "${ENDS_PLUS}"
    awk '$6 == "-"' "${ENDS_FILTERED}"  > "${ENDS_MINUS}"

    # Coverage to bedgraph
    PLUS_BG="${BIGWIG_DIR}/${SRA}.pol.plus.bedgraph"
    MINUS_BG="${BIGWIG_DIR}/${SRA}.pol.minus.bedgraph"

    bedtools genomecov -i "${ENDS_PLUS}"  -g "${CHROM_SIZES}" -bg > "${PLUS_BG}"
    bedtools genomecov -i "${ENDS_MINUS}" -g "${CHROM_SIZES}" -bg > "${MINUS_BG}"

done

