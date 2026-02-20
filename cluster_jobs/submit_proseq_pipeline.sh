#!/usr/bin/env bash
#
# Author: Christos Botos.
# Affiliation: Institute of Molecular Biology and Biotechnology.
# Contact: botoschristos@gmail.com
#
# Script Name: submit_proseq_pipeline.sh
# Description:
#     SLURM job submission script for ProtoPos PRO-seq alignment and
#     processing pipeline on HPC clusters.
#
# Usage:
#     sbatch cluster_jobs/submit_proseq_pipeline.sh
#
# Dependencies:
#     - SLURM workload manager.
#     - conda environment with ProtoPos dependencies (bowtie2, samtools, etc.).
#
# Notes:
#     - Edit the "User-Configurable Section" below before submitting.
#     - SLURM output files are written to cluster_jobs/logs/.
#     - All paths are resolved from this script's location (never hardcoded).
#     - PRO-seq pipelines are I/O-intensive. Ensure sufficient memory for
#       genome indexing and BAM processing.

# ==============================================================================
# SLURM Directives
# ==============================================================================

#SBATCH --job-name=protopos_proseq
#SBATCH --output=cluster_jobs/logs/protopos_proseq_%j.out
#SBATCH --error=cluster_jobs/logs/protopos_proseq_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=CHANGE_ME@example.com

# ==============================================================================
# Path Resolution
# ==============================================================================
# Resolve absolute paths from this script's location (per project policy).

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${script_dir}/.." && pwd)"

# ==============================================================================
# Diagnostics
# ==============================================================================

echo "========================================"
echo "ProtoPos PRO-seq Pipeline — Job Diagnostics"
echo "========================================"
echo "Job ID       : ${SLURM_JOB_ID:-local}"
echo "Job Name     : ${SLURM_JOB_NAME:-local}"
echo "Hostname     : $(hostname)"
echo "Date         : $(date '+%Y-%m-%d %H:%M:%S')"
echo "CPUs         : ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "Memory       : ${SLURM_MEM_PER_NODE:-unknown}"
echo "Working Dir  : ${REPO_DIR}"
echo "========================================"

# ==============================================================================
# Environment Activation
# ==============================================================================
# Auto-detect conda or venv. Edit the environment name if yours differs.

CONDA_ENV_NAME="protopos"

activate_environment() {
    # Try conda first.
    if command -v conda &>/dev/null; then
        echo "Activating conda environment: ${CONDA_ENV_NAME}"
        # Source conda.sh for non-interactive SLURM shells.
        if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
            source "${HOME}/miniconda3/etc/profile.d/conda.sh"
        elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
            source "${HOME}/anaconda3/etc/profile.d/conda.sh"
        elif [[ -f "${HOME}/miniforge3/etc/profile.d/conda.sh" ]]; then
            source "${HOME}/miniforge3/etc/profile.d/conda.sh"
        fi
        conda activate "${CONDA_ENV_NAME}"
        return
    fi

    # Fall back to venv.
    if [[ -d "${REPO_DIR}/venv/bin" ]]; then
        echo "Activating venv: ${REPO_DIR}/venv"
        source "${REPO_DIR}/venv/bin/activate"
        return
    fi

    echo "ERROR: No conda or venv environment found. Exiting."
    exit 1
}

activate_environment

echo "Python       : $(which python)"
echo "Python ver.  : $(python --version)"
echo "========================================"

# ==============================================================================
# User-Configurable Section
# ==============================================================================
# Edit these variables to match your PRO-seq pipeline setup.

PIPELINE_SCRIPT="${REPO_DIR}/scripts/run_pipeline.sh"

# ==============================================================================
# Main Execution
# ==============================================================================

echo ""
echo "Starting ProtoPos PRO-seq pipeline..."
echo "Script       : ${PIPELINE_SCRIPT}"
echo ""

cd "${REPO_DIR}"
bash "${PIPELINE_SCRIPT}"

EXIT_CODE=$?

# ==============================================================================
# Job Summary
# ==============================================================================

echo ""
echo "========================================"
echo "Job Summary"
echo "========================================"
echo "Exit code    : ${EXIT_CODE}"
echo "Completed    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

exit ${EXIT_CODE}
