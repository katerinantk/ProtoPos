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
