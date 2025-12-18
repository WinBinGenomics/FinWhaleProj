#!/bin/bash
#SBATCH --mem=10G
#SBATCH --partition=HIMEM
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# @description  convert a given vcf file to plink format
# Author: Meixi Lin
# Date: Thu Mar 18 14:42:13 2021

# Adjusted by Kaden Winspear -> Loops over correct VCfs to make corresponing plink files.
# usage -> sbatch step3_Prunedvcf2plink.sh

###########################################################
## import packages
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

set -xeo pipefail

WORKDIR="/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/ENP_ESP_admixture"
VCFFILEDIR="/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/ENP_ESP_admixture"
MAFLISTS=('05')  # MAF values to process

# Loop over MAF values
for MAFCUT in "${MAFLISTS[@]}"; do
    # Loop over VCF files but process only the one matching the current MAF
    for VCF in "$VCFFILEDIR"/*.vcf; do
        if [[ -f "$VCF" ]]; then
            # Check if VCF filename contains the current MAF value (e.g., "maf05")
            if [[ "$VCF" == *"maf${MAFCUT}"* ]]; then
                OUT="${WORKDIR}/JointCalls_ESP_ENP_Only_biallelic_all_LDPruned_maf${MAFCUT}_SA_mrf"
                plink --vcf "$VCF" --allow-extra-chr --recode 12 --out "$OUT"
                break  # Exit inner loop after processing the correct VCF
            fi
        else
            echo "Warning: No VCF files found in $VCFFILEDIR."
            exit 1
        fi
    done
done
