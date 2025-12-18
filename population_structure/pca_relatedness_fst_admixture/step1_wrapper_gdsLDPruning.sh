#!/bin/bash
#SBATCH --mem=24G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:00:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL


# @version 		v0
# @usage		sbatch wrapper_step1_gdsLDPruning_all70.sh
# @description	Wrapper to submit the Rscript for LDPruning (NOTE: Using three maf cutoffs)
# Author: Meixi Lin
# Date: Wed Mar 17 11:26:56 2021

# Adjusted / modified by Kaden Winspear @ CSUSM -> Eastern pacific fin whale project.
# Date: April 2025

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
WORKDIR="${HOMEDIR}/fin_genomics/PopStructure/${DATASET}/${REF}"
WORKSCRIPT="${HOMEDIR}/scripts/winfingenomics/population_structure/pca_relatedness_fst_admixture/step1_gdsPASSLDPruning_all70.R"
# The input files
GDSFILE=JointCalls_08_b_VariantFiltration_biallelic_all.gds
OUTPREFIX=JointCalls_all70_filterpass_biallelic_all

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID};"

Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate
