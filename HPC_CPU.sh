#!/bin/bash --login
#SBATCH --job-name=Step00_Check_Resources
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=Step00_%j.log
#SBATCH --error=Step00_%j.err

# Current working directory
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "Working directory: $(pwd)"

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8

echo "Job ID:        ${SLURM_JOB_ID}"
echo "Node:          $(hostname)"
echo "CPUs assigned: ${SLURM_CPUS_PER_TASK}"
echo "Memory:        ${SLURM_MEM_PER_NODE} MB"
echo "Start time:    $(date)"
