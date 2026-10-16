
#!/bin/bash --login
#SBATCH --job-name=Step00_Check_Resources
#SBATCH --account=ascher
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=Step00_%j.log
#SBATCH --error=Step00_%j.err

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8

echo "Working directory: $(pwd)"
echo "Job ID:        ${SLURM_JOB_ID}"
echo "Node:          $(hostname)"
echo "CPUs assigned: ${SLURM_CPUS_PER_TASK}"
echo "Memory:        ${SLURM_MEM_PER_NODE} MB"
echo "Start time:    $(date)"

echo "CPU information:"
lscpu

echo "Memory information:"
free -h

echo "End time: $(date)"
