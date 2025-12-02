#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name=supp_xenium_hlc
#SBATCH --output=slurm_out/supp_xenium_hlc.out
#SBATCH --error=slurm_out/supp_xenium_hlc.err
#SBATCH --time=12:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=12
################################################################################
echo "------------------------------------------------------------------------"
echo "Job Started on $(date)"
echo "System: $(uname -a)"
echo "CPU: $(lscpu | grep 'Model name')"
echo "Total Memory: $(free -h | awk '/Mem:/ {print $2}')"
echo "Node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "------------------------------------------------------------------------"

################################################################################
module load R/4.5.1
module load pandoc/3.2
module load gdal/3.9.0
module load ImageMagick/7.1.1
module load gcc/14.2

export TMPDIR=$(mktemp -d -p /vast/projects/xenium_5k/jazzPanda_paper)
/usr/bin/time -v Rscript -e "rmarkdown::render('supplementary_application_xenium_human_lung_cancer.Rmd', output_format = 'html_document')"
rm -rf  $TMPDIR

echo "------------------------------------------------------------------------"
echo "Job Completed on $(date)"
echo "------------------------------------------------------------------------"
