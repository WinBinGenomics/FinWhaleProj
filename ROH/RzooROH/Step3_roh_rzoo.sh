#! /bin/bash
#SBATCH --job-name=step3_wrapper_roh
#SBATCH --time=90:00:00
#SBATCH --mem=200G
#SBATCH --partition=HIMEM
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=winsp001@csusm.edu

# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Adjusted by Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale project.
# May 1st 2025

# Purpose: This script is a wrapper that runs the scripts for runs of homozygosity.
# Usage: sbatch Step3_roh_rzoo.sh

#Defining directories ---------------------
workdir=/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH

#Main ---------------------------------------
module load R/4.3.1+Bioconductor

set -o pipefail

#Rscript Step3_roh_rzoo_mix10R.R ${workdri}
#Rscript Step3_roh_rzoo_mix15R.R ${workdir}
Rscript Step3_roh_rzoo_mix10R_base3.R ${workdir}
