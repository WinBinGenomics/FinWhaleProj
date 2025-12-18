#!/bin/bash
#SBATCH --account=winsp001
#SBATCH --job-name=pixy_array
#SBATCH --ntasks-per-node=1
#SBATCH --array=01-21
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=22:00:00
#SBATCH --partition=CPU
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

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
VCFDIR="/data/users/snigenda/finwhale_data/blue_ref_omar_data/blue_ref_vcfs"
OUTDIR="/scratch/winsp001/pixy/filtered_VCFs"

mkdir -p ${OUTDIR}

IDX=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})

# Input VCF file
VCF="${VCFDIR}/JointCalls_08_b_VariantFiltration_${IDX}.vcf.gz"

# Output file
OUTVCF="${OUTDIR}/Pixy_filtered_${IDX}.vcf.gz"

# Run filtering
vcftools --gzvcf ${VCF} \
  --remove-indels \
  --max-missing 0.8 \
  --min-meanDP 20 \
  --max-meanDP 500 \
  --recode --stdout | gzip -c > ${OUTVCF}


