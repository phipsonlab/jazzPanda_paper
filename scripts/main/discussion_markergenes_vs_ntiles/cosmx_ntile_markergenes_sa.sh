#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name=ntiles_comp_cosmx_hlc
#SBATCH --output=slurm_out/mg_ntile_cosmx_hlc_%A_%a.out
#SBATCH --error=slurm_out/mg_ntile_cosmx_hlc_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=02:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

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
repeat_times=1
tile_lens="10,20,30,40,50,60,70,80,90,100"

# Pass parameters and SLURM_ARRAY_TASK_ID to the R script
Rscript cosmx_ntile_markergene_result.R $SLURM_ARRAY_TASK_ID $repeat_times $tile_lens


echo "------------------------------------------------------------------------"
echo "Job Completed on $(date)"
echo "------------------------------------------------------------------------"
