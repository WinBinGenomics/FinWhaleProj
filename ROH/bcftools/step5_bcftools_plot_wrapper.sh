#!/bin/bash
#SBATCH --job-name=ROH_plot
#SBATCH --time=24:00:00
#SBATCH --partition=CPU
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu

# Author: Kaden Winspear @ CSUSM. Used to run Rscript for plotting ROH for the populations.

module load R/4.3.1+Bioconductor


script_dir='/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/ROH/bcftools/temp.R'

Rscript ${script_dir}
