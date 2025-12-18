#!/bin/bash
#SBATCH -J Kaden_gone_filter
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=8g
#SBATCH --partition=GPU
#SBATCH --array=01-21
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# Author: Kaden Winspear
# Usage: sbatch step2_FilterVCFs.sh
# This script filters per-population neutral VCFs (already split by pop)
# for: PASS sites, biallelic SNPs, no invariant sites, and no missing data.

module load bcftools/1.19
module load htslib/1.19.1

# Base directory (already split by pop)
base_dir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/split_vcfs"

# Ensure filtered subdirectories exist
mkdir -p "${base_dir}/ENP/filtered" "${base_dir}/ESP/filtered" "${base_dir}/GOC/filtered"

CHR=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})

echo "=== Starting filtering for chromosome ${CHR} ==="

# Filtering function
filter_pop () {
    POP=$1

    IN_VCF="${base_dir}/${POP}/Neutral_sites_SFS_${CHR}.${POP}.vcf.gz"
    OUT_VCF="${base_dir}/${POP}/filtered/Neutral_Pass_BiSnps_${CHR}.${POP}.filtered.vcf.gz"

    echo "[${POP}] Input:  ${IN_VCF}"
    echo "[${POP}] Output: ${OUT_VCF}"

    #bcftools view -f PASS -m2 -M2 -v snps  -Oz -o "${OUT_VCF}" "${IN_VCF}"
    bcftools view  -f PASS -m2 -M2 -v snps -c 1:minor -Oz -o "${OUT_VCF}" "${IN_VCF}"
    tabix -p vcf "${OUT_VCF}"
    echo "[${POP}] Done filtering chromosome ${CHR}"
}

# Run for each population
filter_pop ENP
filter_pop ESP
filter_pop GOC

echo "=== Finished chromosome ${CHR} ==="
