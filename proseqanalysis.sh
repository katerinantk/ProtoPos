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
ALIGN_DIR="alignments"
BT2_INDEX_DIR="bowtie2_index"

GENOME_FA="${DATA_DIR}/genome.fa"
GTF="${DATA_DIR}/gencode.v47.annotation.gtf"

# Bowtie2 index prefix (built from the genome FASTA).
BT2_INDEX="${BT2_INDEX_DIR}/genome"

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
  "${ALIGN_DIR}" \
  "${BT2_INDEX_DIR}"

# -----------------------------
# Build bowtie2 genome index.
# PRO-seq reads nascent unspliced RNA, so a splice-unaware aligner is
# appropriate.
# -----------------------------
if [[ ! -f "${BT2_INDEX}.1.bt2" ]]; then
    echo "Building bowtie2 index..."
    bowtie2-build \
      --threads "${THREADS}" \
      "${GENOME_FA}" \
      "${BT2_INDEX}"
    echo "Bowtie2 index build finished."
else
    echo "Bowtie2 index already exists. Skipping index build."
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
    # Bowtie2 alignment (single-end, unspliced).
    # --very-sensitive for accurate mapping of short nascent RNA fragments.
    # -------------------------
    echo "Aligning ${SRA} with bowtie2..."
    RAW_BAM="${ALIGN_DIR}/${SRA}.raw.bam"
    SORTED_BAM="${ALIGN_DIR}/${SRA}.sorted.bam"

    bowtie2 \
      -x "${BT2_INDEX}" \
      -U "${TRIM}" \
      --very-sensitive \
      --threads "${THREADS}" \
      --no-unal \
    | samtools view -bS -q 10 - \
    | samtools sort -@ "${THREADS}" -o "${SORTED_BAM}" -

    echo "Bowtie2 alignment done for ${SRA}."

    # -------------------------
    # Confirm BAM exists.
    # -------------------------
    if [[ -f "${SORTED_BAM}" ]]; then
        echo "Confirmed: BAM file exists for ${SRA}: ${SORTED_BAM}"
    else
        echo "ERROR: BAM file not found for ${SRA}: ${SORTED_BAM}"
        exit 1
    fi

    # -------------------------
    # Index and deduplicate with umi_tools.
    # PCR duplicates sharing the same UMI and mapping position are collapsed.
    # -------------------------
    DEDUP_BAM="${ALIGN_DIR}/${SRA}.dedup.bam"

    echo "Indexing BAM for ${SRA}..."
    samtools index "${SORTED_BAM}"

    echo "Running UMI deduplication for ${SRA}..."
    umi_tools dedup \
      -I "${SORTED_BAM}" \
      -S "${DEDUP_BAM}" \
      --output-stats="${QC_DIR}/${SRA}.dedup"
    echo "UMI deduplication done for ${SRA}."

    echo "Indexing deduplicated BAM for ${SRA}..."
    samtools index "${DEDUP_BAM}"

done

echo "All samples processed successfully."
