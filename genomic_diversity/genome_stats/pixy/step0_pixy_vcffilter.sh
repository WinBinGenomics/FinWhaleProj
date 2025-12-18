#!/bin/bash
#SBATCH --job-name=pxy_filter
#SBATCH --output=pxyfilter_%j.out
#SBATCH --error=pxyfilter_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=GPU
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


# Author: Kaden Winspear
# Used to filter VCFs prior to running Pixy for simple filtering.
# sbatch [script]

#========WinBin===========
#        .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~
#=======Genomics===========

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate /home/winsp001/.conda/envs/fintools

######## Directories ########
HOMEDIR="/data/shared/snigenda/finwhale_projects/filteredvcf/workshop"
VCFDIR="${HOMEDIR}/vcfs/unfiltered"
OUTDIR="${HOMEDIR}/vcfs/filtered/pixy"

mkdir -p ${OUTDIR}

# In vcf file
VCF="${VCFDIR}/ENP_GOC_chr19-21.vcf.gz"

# out file
outprefix="${OUTDIR}/Passing_biallelic_all_19-21"

# Run filtering
vcftools --gzvcf ${VCF} --remove-filtered-all --remove-indels --max-missing 0.8 --min-meanDP 20 --max-meanDP 400 --recode --out "${outprefix}"

# Compress
bgzip -c "${outprefix}.recode.vcf" > "${outprefix}.vcf.gz"

# Index
tabix -p vcf "${outprefix}.vcf.gz"


