#!/bin/bash
#SBATCH --job-name=rename_plink_pops
#SBATCH --output=rename_plink.%j.out
#SBATCH --error=rename_plink.%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=GPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load mods
module load bcftools/1.19
module load htslib/1.19.1
module load plink/1.9b7.2
module load vcftools/0.1.16


# ***** Dirs + files *****
# Reference dir
REFDIR=/data/shared/snigenda/finwhale_projects/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3

# vcf home dir
VCFDIR=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/split_vcfs

# plink output dirs
PLINKDIR=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/plink
mkdir -p "${PLINKDIR}/ENP" "${PLINKDIR}/ESP" "${PLINKDIR}/GOC"

# autosomes vcfs per population.
ENP_VCF_IN="${VCFDIR}/ENP/filtered/Neutral_Pass_BiSnps_autosomes.ENP.filtered.vcf.gz"
ESP_VCF_IN="${VCFDIR}/ESP/filtered/Neutral_Pass_BiSnps_autosomes.ESP.filtered.vcf.gz"
GOC_VCF_IN="${VCFDIR}/GOC/filtered/Neutral_Pass_BiSnps_autosomes.GOC.filtered.vcf.gz"

# Renamed vcfs output
ENP_VCF_RENAMED="${VCFDIR}/ENP/filtered/Neutral_Pass_BiSnps_autosomes.ENP.filtered.renamedChr.vcf.gz"
ESP_VCF_RENAMED="${VCFDIR}/ESP/filtered/Neutral_Pass_BiSnps_autosomes.ESP.filtered.renamedChr.vcf.gz"
GOC_VCF_RENAMED="${VCFDIR}/GOC/filtered/Neutral_Pass_BiSnps_autosomes.GOC.filtered.renamedChr.vcf.gz"

# ******* Build chromosome mapping from .fai ********

cd "${VCFDIR}"

FAI="${REFDIR}/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta.fai"

if [[ ! -f "$FAI" ]]; then
    echo "ERROR: .fai file not found: $FAI"
    exit 1
fi

echo "Building autosome list from ${FAI} ..."

# Keep only autosomes (exclude sex scaffolds + mt)
grep "NC" "$FAI" | \
  grep -v -e "NC_045806.1" -e "NC_045807.1" -e "NC_001601.1" > autosomesBlueWhale.fai

# autosomes: <scaffold_id> <start> <end>
awk -v OFS="\t" '{print $1, 1, $2}' autosomesBlueWhale.fai > autosomes

# chromosome_lengths
awk -v OFS="\t" '{print FNR, $3}' autosomes > chromosome_lengths.txt

# chromosome_mapping: <old_scaffold_id> <new_chr_number>
# This is what bcftools annotate --rename-chrs uses.
awk -v OFS="\t" '{print $1, FNR}' autosomes > chromosome_mapping

echo "Generated chromosome_mapping (first few lines):"
head chromosome_mapping

# ********** Rename scaffolds 2 chr numbers **********

rename_vcf () {
    local IN_VCF=$1
    local OUT_VCF=$2
    local POP=$3

    echo "----------------------------------------------"
    echo "[${POP}] Renaming chromosomes in VCF:"
    echo "  Input : ${IN_VCF}"
    echo "  Output: ${OUT_VCF}"

    if [[ ! -f "$IN_VCF" ]]; then
        echo "[${POP}] ERROR: Input VCF not found: $IN_VCF"
        exit 1
    fi

    bcftools annotate --rename-chrs chromosome_mapping --output-type z --output "$OUT_VCF" "$IN_VCF"

    if [[ $? -ne 0 ]]; then
        echo "[${POP}] ERROR: bcftools annotate failed."
        exit 1
    fi

    tabix -p vcf "$OUT_VCF"
    echo "[${POP}] Renaming & indexing done."
}

rename_vcf "$ENP_VCF_IN" "$ENP_VCF_RENAMED" "ENP"
rename_vcf "$ESP_VCF_IN" "$ESP_VCF_RENAMED" "ESP"
rename_vcf "$GOC_VCF_IN" "$GOC_VCF_RENAMED" "GOC"

# ******  Convert renamed VCFs to PLINK **********

vcf_to_plink () {
    local IN_VCF=$1
    local POP=$2
    local OUT_PREFIX="${PLINKDIR}/${POP}/${POP}_autosomes_LD"

    echo "----------------------------------------------"
    echo "[${POP}] Converting renamed VCF to PLINK:"
    echo "  VCF: ${IN_VCF}"
    echo "  OUT: ${OUT_PREFIX}.ped / .map"

    vcftools --gzvcf "$IN_VCF" --plink --out "$OUT_PREFIX"

    if [[ $? -ne 0 ]]; then
        echo "[${POP}] ERROR: vcftools --plink failed."
        exit 1
    fi

    # Fix phenotype (6th column) in PED to -9. vcftools does not put this so have to do it manually.
    local PEDFILE="${OUT_PREFIX}.ped"
    local TMP_PED="${OUT_PREFIX}.tmp.ped"

    if [[ ! -f "$PEDFILE" ]]; then
        echo "[${POP}] ERROR: PED file not found after vcftools: $PEDFILE"
        exit 1
    fi

    awk 'BEGIN{OFS=" "} { $6 = -9; print }' "$PEDFILE" > "$TMP_PED" && mv "$TMP_PED" "$PEDFILE"

    echo "[${POP}] PLINK files created and phenotype set to -9."
}

vcf_to_plink "$ENP_VCF_RENAMED" "ENP"
vcf_to_plink "$ESP_VCF_RENAMED" "ESP"
vcf_to_plink "$GOC_VCF_RENAMED" "GOC"

echo "=============================================="
echo "All populations processed: ENP, ESP, GOC."
echo "Renamed VCFs + PLINK files ready for GONE/currentNe."
echo "=============================================="
