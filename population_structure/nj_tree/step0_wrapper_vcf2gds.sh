#!/bin/bash
#SBATCH --mem=40G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# @description	Wrapper to submit the Rscript
# Author: Meixi Lin
# Date: Tue Feb 23 21:42:19 2021

# modified / adjusted by Kaden Winspear @ CSUSM -> Eastern pacific fin whale project. 
# usage: sbatch step0_wrapper_vcf2gds.sh
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
WORKSCRIPT="${HOMEDIR}/scripts/winfingenomics/population_structure/nj_tree/step0_vcf2gds.R"

# The input files
VCFDIR="${HOMEDIR}/filteredvcf/passing_bisnps"
VCFPREFIX="JointCalls_08_b_VariantFiltration"

############################################################
## main

mkdir -p ${WORKDIR}
cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID};"

Rscript --vanilla ${WORKSCRIPT} ${VCFDIR} ${VCFPREFIX} ${WORKDIR}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate
