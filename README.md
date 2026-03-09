# ProtoPos
This repository contains a pipeline for processing PRO-seq nascent transcription data. The workflow downloads raw FASTQ files, performs quality control and trimming, aligns reads with bowtie2 (splice-unaware, appropriate for nascent unspliced RNA), extracts strand-specific 3′ ends and generates genome coverage tracks. It also contains a python pipeline that computes the average PRO-seq signal around the TSS.

## HPC Cluster Submission

Example SLURM job scripts are provided in `cluster_jobs/`:

```bash
# Edit the script to set your email, resource limits, and pipeline script:
vi cluster_jobs/submit_proseq_pipeline.sh

# Submit from the repository root:
sbatch cluster_jobs/submit_proseq_pipeline.sh

# Monitor:
squeue -u $USER
tail -f cluster_jobs/logs/protopos_proseq_<JOB_ID>.out
```

See `cluster_jobs/README.md` for full details.

## CLI Commands

ProtoPos has no pip-installed CLI entry points. The following standalone scripts can be run directly with `python`. All resolve paths relative to their own location via `__file__`, so they work from any directory.

| Script | Description |
|--------|-------------|
| `scripts/gene_histogram.py` | Generate per-gene 3' end distribution histograms from BAM data. |
| `scripts/measure_gene_depth.py` | Measure per-gene sequencing depth from BAM and GTF. |
| `scripts/subset_fastq_by_chrom.py` | Subset FASTQ files by chromosome for test datasets. |
| `scripts/visualize_gene_histograms.py` | Visualize top genes by sequencing depth with histogram plots. |
| `scripts/visualize_gene_histograms_v2.py` | Improved sparse PRO-seq visualizations with smoothing. |
| `scripts/coregulated_profiles.py` | Co-regulated gene profile analysis (22q11.2 region). |
| `scripts/gene_profiles.py` | Random gene profile visualization from BAM data. |
| `scripts/steady_state_figure.py` | Steady-state metagene profile around TSS. |
| `RNAPsignal.py` | RNAP signal extraction around TSS from BigWig data. |

## Coordinate Convention: ProtoPos vs POLYARIS Positions

ProtoPos extracts the **last incorporated nucleotide** (the 3' end of the nascent RNA), which is the position that real sequencing (PRO-seq, GRO-seq) measures. POLYARIS records the **Pol II active center**, which is the next base to be incorporated. These differ by exactly 1 bp:

```
POLYARIS active center:  position P  (0-based gene-relative, half-open)
ProtoPos 3' end:         position P - 1  (last incorporated nucleotide)
```

In genomic coordinates the sign of the offset depends on strand:

| Strand | ProtoPos genomic position | Relationship |
|--------|--------------------------|--------------|
| + | POLYARIS genomic - 1 | ProtoPos is 1 bp upstream |
| - | POLYARIS genomic + 1 | ProtoPos is 1 bp upstream (opposite genomic direction) |

This is expected and biologically correct. When comparing ProtoPos output to POLYARIS `RNA_pol_positions` TSV files, add 1 (gene-relative) to ProtoPos positions to recover the POLYARIS active-center convention.
