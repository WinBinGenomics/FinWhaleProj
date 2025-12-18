#!/usr/bin/env bash
#SBATCH --job-name=subset_esp
#SBATCH --output=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/population_structure/pca_relatedness_fst_admixture/ENP_ESP_admixture/reports/subset_esp.enp.%j.out
#SBATCH --error=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/population_structure/pca_relatedness_fst_admixture/ENP_ESP_admixture/reports/subset_esp.enp%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=HIMEM

module load snigenda-anaconda3/2024.02
source /cm/shared/apps/snigenda-anaconda3/2024.02/etc/profile.d/conda.sh
conda activate GentoolsPop

set -euo pipefail

# ---- user-configurable inputs ----
VCF_IN="/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/all70/Blue/DemoPassInds/DemoPassIds_LDPruned_maf05_SNPs.vcf.gz"
SAMPLE_LIST="/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/population_structure/pca_relatedness_fst_admixture/ENP_ESP_admixture/ENP_ESPids.txt"
VCF_OUT="/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/ENP_ESP_admixture/LDPruned_maf05_SNPs.ESP_ENP_Only.vcf.gz"

# Create logs dir if missing
mkdir -p logs

echo "Job $SLURM_JOB_ID starting on $(hostname)"
echo "VCF_IN: $VCF_IN"
echo "SAMPLE_LIST: $SAMPLE_LIST"
echo "VCF_OUT: $VCF_OUT"

# Ensure input is indexed (force create/overwrite if needed)
bcftools index -f "$VCF_IN"

# Subset to ESP samples
bcftools view -S "$SAMPLE_LIST" -Oz -o "$VCF_OUT" "$VCF_IN"

# Index output
bcftools index -t -f "$VCF_OUT"

echo "Done."
echo "Output: $VCF_OUT"
echo "Index : ${VCF_OUT}.tbi"

