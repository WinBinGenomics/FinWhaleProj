#!/bin/bash
#SBATCH --job-name=ROH_bcftools
#SBATCH --time=24:00:00
#SBATCH --partition=HIMEM
#SBATCH --mem=25G
#SBATCH --array=1-21
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu


# This script identifies regions with autozigosity or runs of homozygosity (ROH) for each population
# Author: Sergio Nigenda & Adjusted by Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale Project.
# Usage: sbatch roh_bcftools.sh population

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

set -o pipefail


########## Defining variables, directories and files

pop=$1

homedir=/data/shared/snigenda/finwhale_projects/fin_genomics
workdir=${homedir}/ROH/rohbcftools
scriptdir=${homedir}/scripts/winfingenomics/ROH/bcftools
vcfdir=${homedir}/ROH/rohbcftools/${pop}
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})
popdir=${workdir}/${pop}

mkdir -p ${workdir}
mkdir -p ${popdir}


########## Performing analysis

cd ${workdir}

bcftools roh -G30 ${vcfdir}/${pop}_PASS_${IDX}.vcf.gz --include 'FILTER="PASS"' -o ${popdir}/${pop}_roh_bcftools_AFACANGT_${IDX}


conda deactivate
