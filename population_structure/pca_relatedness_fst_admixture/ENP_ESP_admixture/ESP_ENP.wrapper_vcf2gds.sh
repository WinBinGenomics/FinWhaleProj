#!/bin/bash
#SBATCH --mem=40G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# Usage for slurm clusters sbatch wrapper_step0_vcf2gds_all50_20210316.sh
# @description	Wrapper to submit the Rscript
# Author: Meixi Lin
# Date: Tue Feb 23 21:42:19 2021

# Adjusted by Kaden Winspear @ CSUSM. Analyzing Eastern pacific fin whales project.
# Date Adjusted: March 5, 2025

############################################################
## import packages

module load R/4.3.1+Bioconductor
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
DATASET='all70'
REF='Blue'

HOMEDIR="/data/shared/snigenda/finwhale_projects"
WORKDIR="${HOMEDIR}/fin_genomics/PopStructure/ENP_ESP_admixture"
WORKSCRIPT="${HOMEDIR}/scripts/winfingenomics/population_structure/pca_relatedness_fst_admixture/ENP_ESP_admixture/ESP_ENP.vcf2gds.R"

####### The input files ########
VCFDIR="${HOMEDIR}/fin_genomics/PopStructure/ENP_ESP_admixture"
VCFPREFIX="LDPruned_maf05_SNPs.ESP_ENP_Only"

############################################################
## main

mkdir -p ${WORKDIR}
cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}."

Rscript --vanilla ${WORKSCRIPT} ${VCFDIR} ${VCFPREFIX} ${WORKDIR}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate
