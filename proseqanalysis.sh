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
  "${STAR_INDEX_DIR}"

# -----------------------------
# Build STAR genome index 
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
# Download / Fastqc / fastp / STAR 
# -----------------------------
for SRA in "${SRA_IDS[@]}"; do
    echo "Processing ${SRA}..."

    RAW="${FASTQ_DIR}/${SRA}.fastq.gz"
    TRIM="${FASTP_DIR}/${SRA}.trimmed.fastq.gz"

    HTML_REPORT="${QC_DIR}/${SRA}.fastp.html"
    JSON_REPORT="${QC_DIR}/${SRA}.fastp.json"

    # -------------------------
    # Download fastq
    # -------------------------
    if [[ ! -f "${RAW}" ]]; then
        echo "Downloading ${SRA}..."
        fasterq-dump "${SRA}" -O "${FASTQ_DIR}" --threads "${THREADS}"
        # if this ever errors, try: fasterq-dump "${SRA}" -O "${FASTQ_DIR}" -e "${THREADS}"
        gzip -f "${FASTQ_DIR}/${SRA}.fastq"
    else
        echo "Skipping download for ${SRA}: ${RAW} already exists."
    fi

    # -------------------------
    # Fastqc on raw reads
    # -------------------------
    echo "Running FastQC on raw reads for ${SRA}..."
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${RAW}"
    echo "FastQC on raw reads done for ${SRA}."

    # -------------------------
    # Trimming with fastp 
    # -------------------------
    echo "Running fastp for ${SRA}..."
    fastp \
      -i "${RAW}" \
      -o "${TRIM}" \
      -w "${THREADS}" \
      --html "${HTML_REPORT}" \
      --json "${JSON_REPORT}"
    echo "fastp done for ${SRA}."

    # -------------------------
    # Fastqc on trimmed reads
    # -------------------------
    echo "Running FastQC on trimmed reads for ${SRA}..."
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${TRIM}"
    echo "FastQC on trimmed reads done for ${SRA}."

    # -------------------------
    # STAR alignment 
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
    # Confirm BAM and convert to SAM
    # -------------------------
    BAM="${STAR_DIR}/${SRA}.Aligned.sortedByCoord.out.bam"
    SAM="${STAR_DIR}/${SRA}.Aligned.sortedByCoord.out.sam"

    if [[ -f "${BAM}" ]]; then
        echo "Confirmed: BAM file exists for ${SRA}: ${BAM}"
    else 
        echo "ERROR: BAM file not found for ${SRA}: ${BAM}"
        exit 1
    fi

    echo "Converting BAM to SAM for ${SRA}..."
    samtools view -h -o "${SAM}" "${BAM}"
    echo "SAM file created for ${SRA}: ${SAM}"

done

echo "All samples processed successfully."
