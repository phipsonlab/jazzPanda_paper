#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name=ntiles_hexbin_comp_sim
#SBATCH --output=slurm_out/ntiles_hexbin_comp_sim_%A_%a.out
#SBATCH --error=slurm_out/ntiles_hexbin_comp_sim_%A_%a.err
#SBATCH --array=1-40
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=6

################################################################################
echo "------------------------------------------------------------------------"
echo "Job Started on $(date)"
echo "------------------------------------------------------------------------"

echo "------------------------------------------------------------------------"
echo "CPU Information:"
lscpu
echo "------------------------------------------------------------------------"

echo "Memory Information:"
free -m
echo "------------------------------------------------------------------------"

echo "System and Kernel:"
uname -a
echo "------------------------------------------------------------------------"
################################################################################
msg "======== begin ========"           #
msg 'Working directory ' `pwd`          # current job working directory
msg 'Job run on nodes ' $SLURM_NODELIST # current job assigned nodes
msg 'Job ntasks assign ' $SLURM_TASKS_PER_NODE #
msg 'NCPU per task' $SLURM_CPUS_PER_TASK # Number of CPUs per task
msg 'Total allocated cores ' $SLURM_CPUS_ON_NODE     # calculated total allocated NCPU
msg 'Job ID ' ${SLURM_JOB_ID}          # Job ID
msg 'Job name ' $SLURM_JOB_NAME        # Job name

################################################################################
module load R/4.5.1
module load pandoc/3.2
module load gdal/3.9.0
module load ImageMagick/7.1.1
module load gcc/14.2

# Define parameters
repeat_times=5
tile_lens="620,310,207,155,124,103,89,78"

export TMPDIR=$(mktemp -d -p .)
# Pass parameters and SLURM_ARRAY_TASK_ID to the R script
Rscript complexity_ntiles_hexbin_simulation.R $SLURM_ARRAY_TASK_ID $repeat_times $tile_lens
rm -rf  $TMPDIR

echo "------------------------------------------------------------------------"
echo "Job Completed on $(date)"
echo "------------------------------------------------------------------------"
