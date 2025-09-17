#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name=run_mg_cosmx_human_healthy_liver_5core
#SBATCH --output=slurm_out/mg_cosmx_human_healthy_liver_5core.out
#SBATCH --error=slurm_out/mg_cosmx_human_healthy_liver_5core.err
#SBATCH --time=12:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=12

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

echo "Disk Space:"
df -h
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

export TMPDIR=$(mktemp -d -p .)
/usr/bin/time -v Rscript mg_cosmx_human_healthy_liver.R
rm -rf  $TMPDIR

echo "------------------------------------------------------------------------"
echo "Job Completed on $(date)"
echo "------------------------------------------------------------------------"
