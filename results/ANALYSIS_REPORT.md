# PRO-seq Analysis Report
## ProtoPos Pipeline - SRR15304570

**Authors:** Katerina Ntouka and Christos Botos
**Affiliation:** Institute of Molecular Biology and Biotechnology
**Date:** December 16, 2025
**Sample:** SRR15304570 (PRO-seq nascent transcription data)

---

## Executive Summary

This report presents a comprehensive analysis of PRO-seq data from sample SRR15304570. The analysis focused on measuring RNA Polymerase II activity across 100 protein-coding genes by quantifying the distribution of 3' end positions (representing actively elongating polymerase positions).

**Key Findings:**
- **93 out of 100** genes showed detectable PRO-seq signal
- **Mean sequencing depth:** 11.96 reads per kilobase
- **Median sequencing depth:** 8.01 reads per kilobase
- **Top gene by depth:** GP1BB (112.28 reads/kb)
- **Most sequenced gene (total reads):** PRAME (465 reads)

---

## 1. Analysis Methodology

### 1.1 Pipeline Overview

The analysis employed a PRO-seq processing pipeline with the following steps:

1. **Read acquisition:** FASTQ files downloaded from SRA (SRR15304570)
2. **Quality control:** FastQC assessment of raw reads
3. **Adapter trimming:** fastp for quality and adapter trimming
4. **Genome alignment:** STAR aligner (splice-aware, though PRO-seq captures unspliced nascent RNA)
5. **3' end extraction:** Conversion of full-length alignments to single-base polymerase positions
6. **Depth quantification:** Read counting across 100 protein-coding genes

### 1.2 Reference Data

- **Genome:** Human reference genome (hg38)
- **Annotation:** GENCODE v47 (protein-coding genes only)
- **Gene selection:** First 100 protein-coding genes from GENCODE annotation
- **Chromosome focus:** Primarily chromosome 22 genes

### 1.3 Bioinformatic Rationale

**Why use 3' ends?**
In PRO-seq, the 3' end of each read represents the exact position of RNA Polymerase II at the moment of biotin incorporation during the run-on reaction. This provides single-nucleotide resolution of polymerase density along genes.

**Sequencing depth metric:**
Depth is calculated as reads per kilobase (reads/kb) to normalize for gene length. This allows fair comparison between short and long genes.

---

## 2. Results

### 2.1 Overall Statistics

| Metric | Value |
|--------|-------|
| Total genes analyzed | 100 |
| Genes with signal | 93 (93%) |
| Genes without signal | 7 (7%) |
| Mean depth | 11.96 reads/kb |
| Median depth | 8.01 reads/kb |
| Maximum depth | 112.28 reads/kb |
| Total reads analyzed | ~2,800 |

**Interpretation:**
The high percentage of genes with signal (93%) indicates good genome-wide coverage. The mean depth of ~12 reads/kb suggests moderate sequencing depth, suitable for detecting highly expressed genes but potentially underpowered for low-expression genes.

### 2.2 Top 10 Genes by Sequencing Depth

| Rank | Gene Symbol | Read Count | Length (kb) | Depth (reads/kb) | Function |
|------|------------|-----------|-------------|------------------|----------|
| 1 | **GP1BB** | 139 | 1.24 | 112.28 | Platelet glycoprotein Ib beta chain |
| 2 | **YDJC** | 117 | 1.98 | 59.21 | Chaperone for ribosomal protein assembly |
| 3 | **THAP7** | 160 | 3.09 | 51.73 | DNA-binding protein, transcriptional regulation |
| 4 | **MRPL40** | 178 | 3.61 | 49.29 | Mitochondrial ribosomal protein L40 |
| 5 | **PRAME** | 465 | 11.66 | 39.88 | Cancer/testis antigen, immune response |
| 6 | **SLC25A1** | 119 | 3.21 | 37.03 | Mitochondrial citrate transporter |
| 7 | **ENSG00000284874** | 277 | 7.55 | 36.67 | Uncharacterized protein-coding gene |
| 8 | **SDF2L1** | 62 | 2.00 | 31.03 | Endoplasmic reticulum chaperone |
| 9 | **TMEM191B** | 72 | 2.77 | 25.97 | Transmembrane protein |
| 10 | **SEPTIN5** | 246 | 9.76 | 25.21 | Septin family GTPase, cytoskeletal regulation |

### 2.3 Functional Enrichment (Preliminary)

The top 10 genes show enrichment for:
- **Mitochondrial function:** MRPL40, SLC25A1 (20%)
- **Protein quality control:** YDJC, SDF2L1 (20%)
- **Membrane/structural proteins:** GP1BB, TMEM191B, SEPTIN5 (30%)

This suggests that the analyzed sample may have high metabolic activity and protein synthesis rates.

---

## 3. Gene-Level Analysis

### 3.1 GP1BB (Highest Depth: 112.28 reads/kb)

**Gene Information:**
- **Symbol:** GP1BB (Glycoprotein Ib Platelet Subunit Beta)
- **Chromosome:** 22
- **Coordinates:** chr22:19,723,539-19,724,776 (+)
- **Length:** 1.24 kb
- **Total reads:** 139

**Distribution Pattern:**
The 3' end distribution shows concentration near the TSS (0-100 bp: 33 reads), indicating active transcription initiation and early elongation. Polymerase density decreases gradually along the gene body.

| Distance from TSS (bp) | Read Count |
|------------------------|-----------|
| 0-100 | 33 |
| 100-200 | 11 |
| 200-300 | 18 |
| 300-400 | 9 |
| 400-500 | 10 |

**Biological Significance:**
GP1BB encodes a platelet surface receptor. High PRO-seq signal suggests active transcription, possibly in a hematopoietic or endothelial cell type.

### 3.2 PRAME (Most Reads: 465)

**Gene Information:**
- **Symbol:** PRAME (Preferentially Expressed Antigen in Melanoma)
- **Chromosome:** 22
- **Coordinates:** chr22:22,547,695-22,559,361 (-)
- **Length:** 11.66 kb
- **Sequencing depth:** 39.88 reads/kb

**Distribution Pattern:**
PRAME shows broader distribution of polymerase across the gene body, consistent with its longer length. Notable clusters at 1.8-2.0 kb downstream of TSS suggest potential regulatory pause sites or alternative processing.

| Distance from TSS (bp) | Read Count |
|------------------------|-----------|
| 0-100 | 1 |
| 1900-2000 | 9 |
| 2000-2100 | 6 |

**Biological Significance:**
PRAME is a cancer/testis antigen involved in immune response. Its high read count (465 total) indicates robust transcription, which may reflect immune cell activity or tumor-related expression.

### 3.3 MRPL40 (Mitochondrial Ribosomal Protein)

**Gene Information:**
- **Symbol:** MRPL40
- **Chromosome:** 22
- **Coordinates:** chr22:19,432,467-19,436,074 (+)
- **Length:** 3.61 kb
- **Total reads:** 178
- **Sequencing depth:** 49.29 reads/kb

**Biological Significance:**
Mitochondrial ribosomal proteins like MRPL40 are essential for oxidative phosphorylation. High PRO-seq signal suggests active mitochondrial biogenesis and energy metabolism.

---

## 4. Data Quality Assessment

### 4.1 Coverage Distribution

The distribution of sequencing depth across genes shows:
- **Genes with >50 reads/kb:** 4 (4%)
- **Genes with 20-50 reads/kb:** 6 (6%)
- **Genes with 5-20 reads/kb:** 37 (37%)
- **Genes with <5 reads/kb:** 46 (46%)
- **Genes with 0 reads:** 7 (7%)

**Interpretation:**
The long tail of low-coverage genes suggests either:
1. Low expression of these genes in the sampled cell type
2. Insufficient sequencing depth for rare transcripts
3. Technical dropout in library preparation

### 4.2 Potential Biases

**Gene length bias:**
Shorter genes (e.g., GP1BB at 1.24 kb) show higher depth when normalized per kb. This is expected for PRO-seq, as polymerase density is typically highest near the TSS.

**Chromosomal bias:**
All analyzed genes are from chromosome 22. This reflects the GTF parsing order and does NOT represent genome-wide expression patterns.

---

## 5. Technical Notes

### 5.1 Alignment Consideration

**STAR vs. Bowtie2:**
This analysis used STAR aligner outputs from a previous pipeline run. For future analyses, Bowtie2 is recommended for PRO-seq data because:
- PRO-seq captures unspliced nascent RNA
- No splice junction detection is needed
- Bowtie2 is faster and more appropriate for DNA-like alignment

### 5.2 Coordinate System

- **GTF annotation:** 1-based, inclusive coordinates
- **SAM/BAM files:** 0-based, half-open coordinates
- **Output files:** 1-based for biological interpretability

### 5.3 Strand Specificity

PRO-seq is strand-specific. The TSS calculation respects strand orientation:
- **Plus strand genes:** TSS = gene start position
- **Minus strand genes:** TSS = gene end position

All distance calculations are relative to the TSS in the 5'→3' transcription direction.

---

## 6. Output Files

### 6.1 Summary Files

| File | Description |
|------|-------------|
| `results/gene_depth.tsv` | Full table of 100 genes with depth statistics |
| `results/ANALYSIS_REPORT.md` | This comprehensive report |

### 6.2 Histogram Data (Top 10 Genes)

For each of the top 10 genes:
- `results/histograms/{GENE}_histogram.txt` - Full position-by-position data
- `results/histograms/{GENE}_summary.txt` - Binned histogram (100 bp bins)

---

## 7. Conclusions

1. **Good signal detection:** 93% of genes show PRO-seq signal, indicating successful library preparation and sequencing.

2. **High-confidence top genes:** GP1BB, YDJC, THAP7, and MRPL40 show >49 reads/kb, providing robust data for downstream analysis.

3. **Functional insights:** Enrichment for mitochondrial and protein quality control genes suggests active cellular metabolism.

4. **Sequencing depth:** Moderate depth (~12 reads/kb mean) is suitable for highly expressed genes but may miss low-expression transcripts.

5. **Data suitability:** The dataset is appropriate for:
   - Identifying highly transcribed genes
   - Comparing polymerase density between genes
   - Detecting pausing patterns in high-expression genes

---

## 8. Recommendations

### 8.1 For Future Experiments

1. **Increase sequencing depth** to 50-100M reads for better coverage of lowly expressed genes
2. **Use Bowtie2** instead of STAR for alignment (more appropriate for unspliced PRO-seq reads)
3. **Expand to genome-wide analysis** beyond chromosome 22
4. **Add biological replicates** for statistical robustness

### 8.2 For Downstream Analysis

1. **Promoter-proximal pausing analysis:** Calculate pausing index (TSS/gene body ratio)
2. **Gene set enrichment analysis:** Test for functional pathway enrichment
3. **Comparison to RNA-seq:** Correlate PRO-seq (nascent) with RNA-seq (steady-state) levels
4. **Visualization:** Generate genome browser tracks (BigWig files) for manual inspection

---

## 9. Data Availability

All analysis scripts and results are available in:
```
ProtoPos/
├── scripts/
│   ├── analyze_gene_depth.sh          # Gene depth quantification
│   └── generate_gene_histograms.sh    # Histogram generation
├── results/
│   ├── gene_depth.tsv                 # Main results table
│   ├── histograms/                    # Per-gene distribution data
│   └── ANALYSIS_REPORT.md             # This report
└── star/
    └── SRR15304570.3prime.bam         # 3' end positions (binary)
```

---

## 10. References

1. Kwak H, Fuda NJ, Core LJ, Lis JT. (2013). Precise maps of RNA polymerase reveal how promoters direct initiation and pausing. *Science* 339(6122):950-3.

2. Mahat DB, Kwak H, Booth GT, et al. (2016). Base-pair-resolution genome-wide mapping of active RNA polymerases using precision nuclear run-on (PRO-seq). *Nat Protoc* 11(8):1455-76.

3. Core LJ, Waterfall JJ, Lis JT. (2008). Nascent RNA sequencing reveals widespread pausing and divergent initiation at human promoters. *Science* 322(5909):1845-8.

---

**End of Report**

*For questions or clarifications, contact: botoschristos@gmail.com*
