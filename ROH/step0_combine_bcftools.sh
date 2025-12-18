#!/bin/bash
#SBATCH --job-name=combine_roh
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu


# @version 		v1
# @usage		sbatch step0_combine_bcftools.sh
# @description	Combine the bcftools output to a shortened version
# Author: Meixi Lin
# adjusted by Kaden Winspear @ CSUSM.
# Date: May 27th

###########################################################
## import packages


module load R/4.5.0+Bioconductor
set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
#HOMEDIR=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/
WORKSCRIPT=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/ROH/step0_combine_bcftools.R
###########################################################

## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}"

Rscript --vanilla ${WORKSCRIPT}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID} Done"
