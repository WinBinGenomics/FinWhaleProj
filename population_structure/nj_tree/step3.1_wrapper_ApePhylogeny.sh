#!/bin/bash
#SBATCH --mem=12G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# @description	Wrapper to submit the Rscript for Ape Trees with bootstraps (NOTE: Using three maf cutoffs) Not plotting anything
# Author: Meixi Lin
# Date: Tue Mar  2 01:57:00 2021
# Usage: sbatch step3.1_wrapper_ApePhylogeny.sh <gdsfile> <outprefix>

# Adjusted by Kaden Winspear -> eastern pacific fin whale project.
# Date: April 2025

############################################################
## import packages

sleep $((RANDOM % 120))

module load R/4.3.1+Bioconductor
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
GDSFILE=${1}
OUTPREFIX=${2}
DATASET='all70'
REF='Blue'
TODAY=$(date "+%Y%m%d")

HOMEDIR="/data/shared/snigenda/finwhale_projects"
WORKDIR="${HOMEDIR}/fin_genomics/PopStructure/${DATASET}/${REF}"
LOGDIR=${WORKDIR}/logs
mkdir -p ${LOGDIR}

WORKSCRIPT="${HOMEDIR}/scripts/winfingenomics/population_structure/nj_tree/step3.1_ApePhylogeny.R"

LOG=${LOGDIR}/step3.1_ApePhylogeny_${OUTPREFIX}_${TODAY}.log

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}; GDSFILE = ${GDSFILE}; OUTPREFIX = ${OUTPREFIX}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}; GDSFILE = ${GDSFILE}; OUTPREFIX = ${OUTPREFIX}" > ${LOG}

Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX} &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}

conda deactivate
