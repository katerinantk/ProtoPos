# Cluster Jobs — ProtoPos

This directory contains SLURM job submission scripts for running ProtoPos PRO-seq alignment pipelines on HPC clusters.

## Contents

| File | Description |
|------|-------------|
| `submit_proseq_pipeline.sh` | PRO-seq alignment and processing pipeline. |
| `logs/` | SLURM output and error files (gitignored). |

## Quick Start

1. **Edit the script** to configure your pipeline:
   ```bash
   vi cluster_jobs/submit_proseq_pipeline.sh
   ```
   - Set `--mail-user` to your email address.
   - Adjust `--cpus-per-task`, `--mem`, and `--time` for your cluster.
   - Set `CONDA_ENV_NAME` if your environment name differs from `protopos`.
   - Set `PIPELINE_SCRIPT` to point to your pipeline script.

2. **Submit the job** from the repository root:
   ```bash
   cd /path/to/ProtoPos
   sbatch cluster_jobs/submit_proseq_pipeline.sh
   ```

3. **Monitor the job:**
   ```bash
   squeue -u $USER
   tail -f cluster_jobs/logs/protopos_proseq_<JOB_ID>.out
   ```

## Default Resources

| Resource | Default | Notes |
|----------|---------|-------|
| CPUs | 8 | Alignment tools (bowtie2, STAR) scale well with threads. |
| Memory | 64G | Genome indexing and BAM processing are memory-intensive. |
| Time | 48:00:00 | Full PRO-seq pipelines with multiple samples take time. |

## Notes

- All paths in the script are resolved relative to the script's location. You can submit from any working directory.
- SLURM `.out` and `.err` files are written to `cluster_jobs/logs/` and are gitignored.
- The script auto-detects conda or venv environments. Source `conda.sh` is handled automatically for non-interactive SLURM shells.
