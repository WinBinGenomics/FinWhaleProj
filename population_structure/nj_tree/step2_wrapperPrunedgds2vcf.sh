#!/bin/bash
#SBATCH --mem=24G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:00:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# @description	Wrapper to submit the Rscript for converting LD pruned gds to vcf files (NOTE: Using three maf cutoffs)
# Author: Meixi Lin
# Date: Thu Mar 18 14:16:53 2021

# Adjusted / modified by Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale project.
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
MAFLISTS=('05' '10' 'NA')
#MAFLISTS=('05')
HOMEDIR="/data/shared/snigenda/finwhale_projects"
WORKDIR="${HOMEDIR}/fin_genomics/PopStructure/${DATASET}/${REF}"
WORKSCRIPT="${HOMEDIR}/scripts/winfingenomics/population_structure/nj_tree/step2_Prunedgds2vcf.R"

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}"

for MAFCUT in ${MAFLISTS[@]}; do
echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Start."
# The input files
OUTPREFIX=JointCalls_filterpass_biallelic_all_LDPruned_maf${MAFCUT}
GDSFILE=${OUTPREFIX}.gds
Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX}
echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Done."
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate
