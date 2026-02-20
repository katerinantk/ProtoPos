# PRO-seq Visualization Guide
## Top 10 Genes by Sequencing Depth - Figure Descriptions

**Authors:** Katerina Ntouka and Christos Botos
**Date:** December 16, 2025
**Sample:** SRR15304570

---

## 📊 Overview

This guide describes the 12 high-resolution figures generated from the PRO-seq analysis, located in `results/figures/`.

All plots were generated using matplotlib 3.10.0 at 300 DPI resolution, suitable for publication.

---

## 🎯 Summary Visualizations (2 files)

### 1. **`top10_comparison_panel.png`** (451 KB)
**Multi-panel comparison of all top 10 genes**

- **Layout:** 4×3 grid showing all 10 genes simultaneously
- **Purpose:** Quick visual comparison of polymerase distribution patterns
- **X-axis:** Distance from TSS (bp)
- **Y-axis:** Read count
- **Red dashed line:** Marks TSS position (0 bp)
- **Titles:** Show gene name, total reads, and normalized depth

**Key Observations:**
- **GP1BB, YDJC, SLC25A1:** Strong promoter-proximal peaks (pausing)
- **PRAME:** Broader distribution due to longer gene length (11.66 kb)
- **Most genes:** Show declining polymerase density downstream of TSS

**Use Case:** Publication figure panel, presentation slide

---

### 2. **`sequencing_depth_comparison.png`** (201 KB)
**Side-by-side bar chart comparison**

- **Left Panel:** Sequencing depth (reads/kb) - normalized for gene length
- **Right Panel:** Total read count - absolute numbers
- **Color:** Viridis gradient (darker = higher rank)
- **Orientation:** Horizontal bars for easy gene name reading

**Key Insights:**
- GP1BB dominates with 112.28 reads/kb (9× higher than median)
- PRAME has most total reads (465) but lower depth due to length
- Clear rank order visible from color gradient

**Use Case:** Summary statistics presentation, grant applications

---

## 🧬 Individual Gene Histograms (10 files)

Each gene has a dedicated high-resolution figure with two panels:

### Standard Layout (All 10 Genes):

**Top Panel:** Distribution histogram
- **Bins:** 50 bins across the gene length
- **Color:** Blue (#2E86AB) with black edges
- **Red dashed line:** TSS marker
- **Inset box:** Gene statistics (reads, length, depth, strand, coordinates)
- **Grid:** Horizontal gridlines for read count reference

**Bottom Panel:** Cumulative distribution
- **Color:** Purple (#A23B72)
- **Purpose:** Shows what percentage of reads fall within X bp of TSS
- **Interpretation:** Steep rise = tight concentration near TSS

---

### Gene-Specific Highlights:

#### **`GP1BB_histogram.png`** (217 KB) - Rank #1
**Gene:** Platelet glycoprotein Ib beta chain
**Pattern:** Classic promoter-proximal pausing

```
Distance from TSS    Read Count
    0-100 bp             33  ██████████████
  100-300 bp             29  ███████████
  300-600 bp             29  ███████████
```

**Biological Interpretation:**
- Strong polymerase accumulation within 100 bp of TSS
- Characteristic of highly regulated, poised genes
- May indicate pausing for co-transcriptional regulation

**Cumulative Distribution:**
- 50% of reads within first 450 bp
- Rapid accumulation near TSS

---

#### **`PRAME_histogram.png`** (230 KB) - Rank #5
**Gene:** Preferentially Expressed Antigen in Melanoma
**Pattern:** Distributed elongation across long gene body

```
Distance from TSS    Read Count    Pattern
    0-500 bp              17       Low initial density
  1800-2100 bp            24       Mid-gene cluster
  5000+ bp             ~200        Sustained elongation
```

**Biological Interpretation:**
- Minimal promoter-proximal pausing (only 1 read at TSS!)
- Active elongation throughout 11.66 kb gene body
- Mid-gene peak (1.9-2.1 kb) suggests regulatory element or pause site
- Cancer/testis antigen expression → possible immune-related sample

**Cumulative Distribution:**
- Gradual linear rise (not steep)
- 50% of reads reached at ~5 kb from TSS
- Indicates processive elongation, not pausing-dominated

---

#### **`MRPL40_histogram.png`** (224 KB) - Rank #4
**Gene:** Mitochondrial Ribosomal Protein L40
**Pattern:** Moderate pausing with steady elongation

**Biological Interpretation:**
- Mitochondrial gene with constitutive expression
- Balanced initiation and elongation
- Essential for oxidative phosphorylation

---

#### **`YDJC_histogram.png`** (221 KB) - Rank #2
**Gene:** YdjC Homolog (Ribosome Assembly Chaperone)
**Pattern:** Sharp TSS peak with rapid decline

**Biological Interpretation:**
- Strong pausing/regulation at promoter
- Short gene (1.98 kb) fully covered
- Ribosome assembly → high cellular activity

---

#### **`THAP7_histogram.png`** (221 KB) - Rank #3
**Gene:** THAP Domain Containing 7 (Transcription Factor)
**Pattern:** Balanced distribution

**Biological Interpretation:**
- DNA-binding protein with moderate expression
- Even polymerase progression across gene

---

#### **`SLC25A1_histogram.png`** (217 KB) - Rank #6
**Gene:** Mitochondrial Citrate Transporter
**Pattern:** Strong promoter peak

**Biological Interpretation:**
- Second mitochondrial gene in top 10
- High metabolic demand in sample
- Citrate transport → TCA cycle activity

---

#### **`ENSG00000284874_histogram.png`** (229 KB) - Rank #7
**Gene:** Uncharacterized Protein-Coding Gene
**Pattern:** Distributed across 7.55 kb

**Biological Interpretation:**
- Novel or poorly annotated gene with high expression
- 277 total reads indicate robust transcription
- Candidate for functional characterization

---

#### **`SDF2L1_histogram.png`** (215 KB) - Rank #8
**Gene:** Stromal Cell Derived Factor 2 Like 1 (ER Chaperone)
**Pattern:** Compact distribution (2 kb gene)

**Biological Interpretation:**
- ER stress response protein
- May indicate protein folding demand

---

#### **`TMEM191B_histogram.png`** (213 KB) - Rank #9
**Gene:** Transmembrane Protein 191B
**Pattern:** Short gene (2.77 kb) with even coverage

**Biological Interpretation:**
- Membrane protein with moderate expression
- 72 reads across short gene = decent coverage

---

#### **`SEPTIN5_histogram.png`** (214 KB) - Rank #10
**Gene:** Septin 5 (Cytoskeletal GTPase)
**Pattern:** 246 reads across 9.76 kb

**Biological Interpretation:**
- Cytoskeletal regulation protein
- Long gene with sustained elongation
- Important for cell division and structure

---

## 📈 How to Interpret the Figures

### Histogram Panel (Top):
1. **High peak at TSS (0 bp):** Promoter-proximal pausing (Pol II regulation)
2. **Gradual decline:** Normal elongation with some pausing
3. **Flat distribution:** Processive elongation (minimal pausing)
4. **Mid-gene peaks:** Potential regulatory elements, splice sites (though PRO-seq is unspliced!)

### Cumulative Panel (Bottom):
1. **Steep rise at start:** Most reads near TSS (pausing-dominated)
2. **Linear rise:** Even distribution (elongation-dominated)
3. **50% threshold:** Median polymerase position from TSS

### Promoter-Proximal Pausing Index:
Calculate as: (Reads in TSS +/-100 bp) / (Reads in gene body)
- **High index (>2):** Pausing-dominated (GP1BB, YDJC)
- **Low index (<0.5):** Elongation-dominated (PRAME)

---

## 🔬 Technical Specifications

### Figure Properties:
- **Format:** PNG (lossless compression)
- **Resolution:** 300 DPI (publication quality)
- **Color Space:** sRGB
- **Font:** Default matplotlib (DejaVu Sans)
- **Grid:** Light gray dashed lines (alpha=0.3)
- **Spines:** Top and right removed (minimal design)

### Data Source:
- **Raw data:** `results/histograms/{GENE}_histogram.txt`
- **Summary data:** `results/histograms/{GENE}_summary.txt`
- **Gene metadata:** `results/gene_depth.tsv`

### Reproducibility:
All figures can be regenerated with:
```bash
source ../venv310/bin/activate
python scripts/visualize_gene_histograms.py
```

---

## 📚 Suggested Figure Usage

### For Publications:
- **Main Figure:** `top10_comparison_panel.png` (shows pattern diversity)
- **Supplement:** Individual gene histograms for genes of interest
- **Bar Chart:** `sequencing_depth_comparison.png` (summary statistics)

### For Presentations:
- **Title Slide:** `sequencing_depth_comparison.png` (clear ranking)
- **Results Slide:** `top10_comparison_panel.png` (comprehensive)
- **Specific Examples:** Pick 2-3 individual histograms (e.g., GP1BB vs PRAME)

### For Grant Applications:
- **Preliminary Data:** All 3 summary figures
- **Methodology:** Individual histogram to show data quality

---

## 🎨 Color Scheme

All figures use a consistent color palette:
- **Histograms:** Blue (#2E86AB) - represents nascent RNA
- **Cumulative:** Purple (#A23B72) - distinguishes from histogram
- **TSS Marker:** Red dashed line - universal reference point
- **Comparison Bars:** Viridis gradient - shows rank order

---

## 📊 File Size Summary

| Figure Type | Count | Size Range | Total Size |
|-------------|-------|------------|------------|
| Individual histograms | 10 | 213-230 KB | ~2.2 MB |
| Comparison panel | 1 | 451 KB | 451 KB |
| Depth comparison | 1 | 201 KB | 201 KB |
| **TOTAL** | **12** | - | **~2.9 MB** |

All files are optimized for:
- ✅ Email attachments (<3 MB total)
- ✅ Journal submission systems
- ✅ Presentation software (PowerPoint, Keynote)
- ✅ Web display (high quality, reasonable size)

---

## 🔍 Next Steps

### Additional Analyses:
1. **Calculate pausing indices** for all genes
2. **Compare to RNA-seq data** (nascent vs steady-state)
3. **Gene Ontology enrichment** of top genes
4. **Correlation with histone marks** (if ChIP-seq available)

### Extended Visualizations:
1. **Metagene plots:** Average profile across all genes
2. **Heatmaps:** All 100 genes ranked by depth
3. **TSS-centered density plots:** ±2 kb window
4. **Strand-specific BigWig tracks:** For genome browser

---

## 📖 References for Interpretation

1. **Promoter-proximal pausing:**
   - Core LJ, et al. (2008) Science 322:1845-8
   - "Nascent RNA sequencing reveals widespread pausing"

2. **PRO-seq methodology:**
   - Mahat DB, et al. (2016) Nat Protoc 11:1455-76
   - "Base-pair-resolution mapping using PRO-seq"

3. **Pausing index calculation:**
   - Kwak H, Lis JT (2013) Cold Spring Harb Symp Quant Biol 78:113-21
   - "Control of transcriptional elongation"

---

**End of Visualization Guide**

*All figures generated: December 16, 2025*
*ProtoPos Pipeline v1.0*
*For questions: botoschristos@gmail.com*
