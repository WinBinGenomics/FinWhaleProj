#!/bin/bash
#SBATCH -J Kaden_concat_vcfs
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=8g
#SBATCH --partition=GPU
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# Author: Kaden Winspear
# Usage: sbatch step3_ConcatVCFs.sh
# Goal: Concatenate per-chromosome filtered VCFs into one whole-genome VCF per population.
# Result: 3 VCFs (ENP, ESP, GOC) + indices.

module load bcftools/1.19
module load htslib/1.19.1

base_dir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/split_vcfs"

# Output files (whole-genome per pop)
OUT_ENP="${base_dir}/ENP/filtered/Neutral_Pass_BiSnps_autosomes.ENP.filtered.vcf.gz"
OUT_ESP="${base_dir}/ESP/filtered/Neutral_Pass_BiSnps_autosomes.ESP.filtered.vcf.gz"
OUT_GOC="${base_dir}/GOC/filtered/Neutral_Pass_BiSnps_autosomes.GOC.filtered.vcf.gz"

concat_pop () {
    local POP=$1
    local OUT_VCF=$2

    echo "=== Concatenating population: ${POP} ==="

    local files=()
    for CHR in $(seq -w 1 21); do
        local f="${base_dir}/${POP}/filtered/Neutral_Pass_BiSnps_${CHR}.${POP}.filtered.vcf.gz"
        if [[ ! -f "$f" ]]; then
            echo "[ERROR] Missing file: $f"
            exit 1
        fi
        files+=("$f")
    done

    echo "[${POP}] Files to concat:"
    printf '  %s\n' "${files[@]}"

    # Concatenate in chrom order
    bcftools concat  -Oz -o "${OUT_VCF}" "${files[@]}"

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] bcftools concat failed for ${POP}"
        exit 1
    fi

    tabix -p vcf "${OUT_VCF}"
    echo "[${POP}] Done. Output: ${OUT_VCF}"
}

concat_pop ENP "${OUT_ENP}"
concat_pop ESP "${OUT_ESP}"
concat_pop GOC "${OUT_GOC}"

echo "=== All populations concatenated successfully ==="
