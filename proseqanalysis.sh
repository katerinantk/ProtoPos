#!/usr/bin/env bash

# -----------------------------
# Inputs
# -----------------------------
SRA_IDS=("SRR15304569" "SRR15304570")
THREADS=8

# UMI pattern at the 5' end of R1 (6 random bases in the standard Mahat et al.
# 2016 PRO-seq protocol).  Verify this matches the protocol used for your data.
UMI_PATTERN="NNNNNN"

# Standard Illumina TruSeq small RNA 3' adapter used in most PRO-seq protocols.
# Verify this matches the adapter used for your data.
ADAPTER_SEQ="TGGAATTCTCGGGTGCCAAGG"

DATA_DIR="data"
FASTQ_DIR="fastq"
UMI_DIR="umi"
FASTP_DIR="fastp"
QC_DIR="qc"
STAR_DIR="star"
STAR_INDEX_DIR="starindex"

GENOME_FA="${DATA_DIR}/genome.fa"
GTF="${DATA_DIR}/gencode.v47.annotation.gtf"

# -----------------------------
# Check if the FASTA and GTF files exist.
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
# Create output directories.
# -----------------------------
mkdir -p \
  "${FASTQ_DIR}" \
  "${UMI_DIR}" \
  "${FASTP_DIR}" \
  "${QC_DIR}" \
  "${STAR_DIR}" \
  "${STAR_INDEX_DIR}"

# -----------------------------
# Build STAR genome index.
# -----------------------------
if [[ ! -f "${STAR_INDEX_DIR}/Genome" ]]; then
    echo "Building STAR index..."
    STAR \
      --runThreadN "${THREADS}" \
      --runMode genomeGenerate \
      --genomeDir "${STAR_INDEX_DIR}" \
      --genomeFastaFiles "${GENOME_FA}" \
      --sjdbGTFfile "${GTF}" \
      --sjdbOverhang 75
    echo "STAR index build finished."
else
    echo "STAR index already exists. Skipping index build."
fi

# -----------------------------
# Download / UMI extract / Trim / Align / Dedup
# -----------------------------
for SRA in "${SRA_IDS[@]}"; do
    echo "Processing ${SRA}..."

    RAW="${FASTQ_DIR}/${SRA}.fastq.gz"
    UMI_EXTRACTED="${UMI_DIR}/${SRA}.umi.fastq.gz"
    TRIM="${FASTP_DIR}/${SRA}.trimmed.fastq.gz"

    HTML_REPORT="${QC_DIR}/${SRA}.fastp.html"
    JSON_REPORT="${QC_DIR}/${SRA}.fastp.json"

    # -------------------------
    # Download fastq.
    # -------------------------
    if [[ ! -f "${RAW}" ]]; then
        echo "Downloading ${SRA}..."
        fasterq-dump "${SRA}" -O "${FASTQ_DIR}" --threads "${THREADS}"
        # If this ever errors, try: fasterq-dump "${SRA}" -O "${FASTQ_DIR}" -e "${THREADS}"
        gzip -f "${FASTQ_DIR}/${SRA}.fastq"
    else
        echo "Skipping download for ${SRA}: ${RAW} already exists."
    fi

    # -------------------------
    # FastQC on raw reads.
    # -------------------------
    echo "Running FastQC on raw reads for ${SRA}..."
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${RAW}"
    echo "FastQC on raw reads done for ${SRA}."

    # -------------------------
    # Extract UMIs from the 5' end of reads.
    # The UMI is moved from the read sequence into the read name so that
    # umi_tools dedup can use it after alignment.
    # -------------------------
    echo "Extracting UMIs for ${SRA} (pattern: ${UMI_PATTERN})..."
    umi_tools extract \
      --stdin="${RAW}" \
      --stdout="${UMI_EXTRACTED}" \
      --bc-pattern="${UMI_PATTERN}"
    echo "UMI extraction done for ${SRA}."

    # -------------------------
    # Trimming with fastp.
    # Explicit adapter sequence for the standard PRO-seq 3' RNA adapter.
    # -------------------------
    echo "Running fastp for ${SRA}..."
    fastp \
      -i "${UMI_EXTRACTED}" \
      -o "${TRIM}" \
      -w "${THREADS}" \
      --adapter_sequence "${ADAPTER_SEQ}" \
      --html "${HTML_REPORT}" \
      --json "${JSON_REPORT}"
    echo "fastp done for ${SRA}."

    # -------------------------
    # FastQC on trimmed reads.
    # -------------------------
    echo "Running FastQC on trimmed reads for ${SRA}..."
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${TRIM}"
    echo "FastQC on trimmed reads done for ${SRA}."

    # -------------------------
    # STAR alignment.
    # -------------------------
    echo "Aligning ${SRA} with STAR..."
    PREFIX="${STAR_DIR}/${SRA}."
    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX_DIR}" \
      --readFilesIn "${TRIM}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${PREFIX}" \
      --outSAMtype BAM SortedByCoordinate
    echo "STAR alignment done for ${SRA}."

    # -------------------------
    # Confirm BAM exists.
    # -------------------------
    BAM="${STAR_DIR}/${SRA}.Aligned.sortedByCoord.out.bam"

    if [[ -f "${BAM}" ]]; then
        echo "Confirmed: BAM file exists for ${SRA}: ${BAM}"
    else
        echo "ERROR: BAM file not found for ${SRA}: ${BAM}"
        exit 1
    fi

    # -------------------------
    # Index and deduplicate with umi_tools.
    # PCR duplicates sharing the same UMI and mapping position are collapsed.
    # -------------------------
    DEDUP_BAM="${STAR_DIR}/${SRA}.Aligned.sortedByCoord.out.dedup.bam"

    echo "Indexing BAM for ${SRA}..."
    samtools index "${BAM}"

    echo "Running UMI deduplication for ${SRA}..."
    umi_tools dedup \
      -I "${BAM}" \
      -S "${DEDUP_BAM}" \
      --output-stats="${QC_DIR}/${SRA}.dedup"
    echo "UMI deduplication done for ${SRA}."

    echo "Indexing deduplicated BAM for ${SRA}..."
    samtools index "${DEDUP_BAM}"

done

echo "All samples processed successfully."
