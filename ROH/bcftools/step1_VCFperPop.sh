#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --mem=15G
#SBATCH --job-name=NewVcf_bcftools
#SBATCH --array=1-21
#SBATCH --mail-type=ALL

# This script create new vcf files per population that later will be used as input in runs of homozygosity (ROH) analysis
# Author: Sergio Nigenda (UCLA) & Kaden Winspear @ CSUSM.
# June 2nd 2025
# Usage: sbatch newVcf_perPop.sh population


module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

set -o pipefail


########## Defining variables, directories and files

pop=$1

homedir=/data/shared/snigenda/finwhale_projects/fin_genomics
workdir=${homedir}/ROH/rohbcftools/
scriptdir=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/ROH/bcftools
vcfdir=/data/shared/snigenda/finwhale_projects/filteredvcf/passing_bisnps/
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})
popdir=${workdir}/${pop}
popfile=${scriptdir}/${pop}_samples.txt

mkdir -p ${workdir}
mkdir -p ${popdir}


########## Performing analysis

cd ${workdir}

# The bcftools manual recomends to do the subsetting and filtering of as two different steps becuase sample removal can change the allele frequency, so we do it like it is recommended by pipping the results of step 1 (subsetting samples) into step 2 (filtering variants) 
bcftools view -S ${popfile} \
-Ou ${vcfdir}/JointCalls_08_b_VariantFiltration_${IDX}_passing_bisnps.vcf.gz | bcftools view \
--include 'FILTER=="PASS"' \
-o ${popdir}/${pop}_PASS_${IDX}.vcf.gz \
-Oz


conda deactivate
