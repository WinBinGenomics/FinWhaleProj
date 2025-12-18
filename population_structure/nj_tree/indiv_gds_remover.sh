#!/bin/bash
#SBATCH --job-name=PopStructure      # Job name
#SBATCH --mem=20G                     # Total memory request (matches h_data=20G)
#SBATCH --time=23:00:00               # Wall time limit (HH:MM:SS)
#SBATCH --mail-type=ALL    # Notification events (abe = BEGIN,END,FAIL)
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --partition=CPU

# Author: Kaden Winspear @ CSUSM -> Eastern pacific fin whale project.
# Date: April 2024
# Usage: sbatch indiv_gds_remover.sh
# Purpose: uses GATK and SNP relate to remove certain individuals. Used to remove individuals of low genotype.

module load bcftools/1.19
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate whalinGATK

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
DATASET='all70'
REF='Blue'
MAFCUT='05'
TODAY=$(date "+%Y%m%d")

HOMEDIR=/data/shared/snigenda/finwhale_projects/fin_genomics
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}
OUTDIR=${WORKDIR}/DemoPassInds
mkdir -p ${OUTDIR}

# Blue whale reference genome
REFERENCE=/data/shared/snigenda/finwhale_projects/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta

# samples to exclude
#EXCLUDE_SAMPLE=(ENPOR12" "GOC077") # For LowGenoRm dataset.
EXCLUDE_SAMPLE=("ENPOR12" "ENPCA01" "ENPCA09" "GOC010" "GOC080" "GOC111" "GOC077" "GOC053") # Removed Individuals for Demographic analysis.

# Original VCF file
INVCF="${WORKDIR}/all70Inds/JointCalls_all70_filterpass_biallelic_all_LDPruned_maf05_SA_mrF.vcf"
# Out prefix
OUTPREFIX="DemoPassIds_LDPruned_maf05_SNPs"

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}; Maf cutoff = ${MAFCUT}; Input VCF = ${INVCF}"

############################################################
# Rewrites VCF files without the low genotype Individuals.
cd $OUTDIR

gatk --java-options "-Xmx4G" SelectVariants \
  -R $REFERENCE \
  $(for ss in "${EXCLUDE_SAMPLE[@]}"; do echo "--exclude-sample-name $ss"; done) \
  --exclude-filtered \
  --remove-unused-alternates \
  --exclude-non-variants \
  --restrict-alleles-to BIALLELIC \
  --select-type-to-include SNP \
  -V "${INVCF}" \
  -O "${OUTPREFIX}.vcf"

bcftools +counts ${OUTPREFIX}.vcf

conda deactivate

source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

# Compress and index the VCF for efficient processing
bgzip ${OUTPREFIX}.vcf
tabix ${OUTPREFIX}.vcf.gz

conda deactivate

# Load R module and convert VCF to GDS

module load R/4.3.1+Bioconductor

# Convert VCF to GDS using SNPRelate
Rscript --vanilla -e "
library(gdsfmt)
library(SNPRelate)
vcf_file <- '${OUTDIR}/${OUTPREFIX}.vcf.gz'
gds_file <- '${OUTDIR}/${OUTPREFIX}.gds'
snpgdsVCF2GDS(vcf_file, gds_file, method='biallelic.only')
snpgdsSummary(gds_file)
"
