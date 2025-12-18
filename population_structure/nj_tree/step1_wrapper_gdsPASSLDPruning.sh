#!/bin/bash
#SBATCH --mem=24G
#SBATCH --partition=HIMEM
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# @description	Wrapper to submit the Rscript for LDPruning (NOTE: Using three maf cutoffs)
# Author: Meixi Lin
# Date: Thu Feb 25 00:20:23 2021

# Adjusted by Kaden Winspear (csusm) -> eastern pacific fin whale project.
# sbatch step1_wrapper_gdsPASSLDPruning.sh

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

HOMEDIR='/data/shared/snigenda/finwhale_projects'
WORKDIR="${HOMEDIR}/fin_genomics/PopStructure/${DATASET}/${REF}"
WORKSCRIPT="${HOMEDIR}/scripts/winfingenomics/population_structure/nj_tree/step1_gdsPASSLDPruning.R"

# The input files
GDSFILE=JointCalls_08_b_VariantFiltration_biallelic_all.gds
OUTPREFIX=JointCalls_filterpass_biallelic_all

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID $SLURM_{JOB_ID};"

Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate
