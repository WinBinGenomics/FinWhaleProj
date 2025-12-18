#!/bin/bash
#SBATCH -J Kaden_gone
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=8g
#SBATCH --partition=GPU
#SBATCH --array=01-21
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# Author: Kaden Winspear
# Usage: sbatch step1_SplitPops.sh [Thhis array parallelizes the vcfs per chrom (1-21 in this case)]
# This script is the first step in the LD analysis pipeline. Before any population filtering we first seperate the populations. Admixed and poor genotype IDs are not included.

# Load the mods
module load bcftools/1.19
module load htslib/1.19.1

# inputdir
input_dir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/Neutral_regions/neutralVCFs"

# outputdir
output_dir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/split_vcfs"


# Sample lists for each pop. Admixed and lowgenotype excluded.
ENP_SAMPLES=(ENPAK19 ENPAK20 ENPAK21 ENPAK22 ENPAK23 ENPAK24 ENPAK25 ENPAK26 ENPAK27 ENPAK28 ENPAK29 ENPAK30
  ENPBC16 ENPBC17 ENPBC18 ENPCA02 ENPCA03 ENPCA04 ENPCA05 ENPCA06 ENPCA07 ENPCA08 ENPOR10 ENPOR11 
  ENPOR13 ENPWA14 ENPWA15)

ESP_SAMPLES=(ESPCL01 ESPCL02 ESPCL03 ESPCL04 ESPCL05 ESPCL06 ESPCL07 ESPCL08 ESPCL09 ESPCL10 ESPCL11 ESPCL12
  ESPCL13 ESPCL14 ESPCL15 ESPCL16 ESPCL17 ESPCL18 ESPCL19 ESPCL20)

GOC_SAMPLES=(GOC002 GOC006 GOC025 GOC038 GOC050 GOC063 GOC068 GOC071 GOC082
  GOC086 GOC091 GOC100 GOC112 GOC116 GOC125)

#(GOC002, GOC025, GOC038, GOC050, GOC068, GOC082, GOC086, GOC091, GOC100, GOC112)

# Convert sample lists to comma-separated strings
ENP_LIST=$(IFS=,; echo "${ENP_SAMPLES[*]}")
ESP_LIST=$(IFS=,; echo "${ESP_SAMPLES[*]}")
GOC_LIST=$(IFS=,; echo "${GOC_SAMPLES[*]}")

# Get the VCF file for this task.
VCF_FILE="${input_dir}/Neutral_sites_SFS_$(printf "%02d" $SLURM_ARRAY_TASK_ID).vcf.gz"
# Define output file prefixes.
OUT_PREFIX_ENP="${output_dir}/ENP/Neutral_sites_SFS_$(printf "%02d" $SLURM_ARRAY_TASK_ID).ENP.vcf.gz"
OUT_PREFIX_ESP="${output_dir}/ESP/Neutral_sites_SFS_$(printf "%02d" $SLURM_ARRAY_TASK_ID).ESP.vcf.gz"
OUT_PREFIX_GOC="${output_dir}/GOC/Neutral_sites_SFS_$(printf "%02d" $SLURM_ARRAY_TASK_ID).GOC.vcf.gz"

# Filter vcf by pop
echo "Filtering ENP samples from $VCF_FILE..."
bcftools view -Oz -s "$ENP_LIST" -o "$OUT_PREFIX_ENP" "$VCF_FILE"

echo "Filtering ESP samples from $VCF_FILE..."
bcftools view -Oz -s "$ESP_LIST" -o "$OUT_PREFIX_ESP" "$VCF_FILE"

echo "Filtering GOC samples from $VCF_FILE..."
bcftools view -Oz -s "$GOC_LIST" -o "$OUT_PREFIX_GOC" "$VCF_FILE"

# Index the output VCFs
echo "Indexing output VCFs..."
tabix -p vcf "$OUT_PREFIX_ENP"
tabix -p vcf "$OUT_PREFIX_ESP"
tabix -p vcf "$OUT_PREFIX_GOC"

#EXCLUDE_SAMPLE=(ENPOR12" "ENPCA01" "ENPCA09" "GOC010" "GOC080" "GOC111" "GOC077" "GOC053")
