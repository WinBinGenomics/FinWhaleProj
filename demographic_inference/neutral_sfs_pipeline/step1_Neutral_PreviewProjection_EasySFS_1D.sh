#!/bin/bash
#SBATCH --job-name=PreviewProjection
#SBATCH --time=20:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-21
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu

# Usage @ CSUSM: sbatch easySFS_1_ProjectionPreview_onlyPASSpopHet75.sh
# This script will get and parse project preview for a given vcf file
# @modification Mon Jul 13 21:44:56 2020
# @modification Final settings. Add per population maxHet filters at 0.75.

set -eo pipefail

###########################################################
## import packages

module load anaconda3/2023.09
module load R/4.5.0+Bioconductor
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate whalinGATK

###########################################################
## define functions

# select variants
# $1 = input VCF name
# $2 = output VCF name
gatk_select_variants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start gatk4 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, select SNPs for ${1}" >> ${PROGRESSLOG}
gatk --java-options "-Xmx4G" SelectVariants \
-R $REFERENCE \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude-sample-name ${ss} "; done) \
--exclude-filtered \
--remove-unused-alternates \
--exclude-non-variants \
--restrict-alleles-to BIALLELIC \
--select-type-to-include SNP \
-V ${1} \
-O ${2} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Finish gatk4 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, select SNPs for ${1}" >> ${PROGRESSLOG}
}

# select invariant sites
# $1 = input VCF name
# $2 = output VCF name
gatk_select_invariants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start GATK4 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, Select Invariant sites for ${1}" >> ${PROGRESSLOG}
gatk --java-options "-Xmx4G" SelectVariants \
-R $REFERENCE \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude-sample-name ${ss} "; done) \
--exclude-filtered \
--remove-unused-alternates \
--select-type-to-include NO_VARIATION \
--select-type-to-exclude INDEL \
-V ${1} \
-O ${2} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Finish GATK4 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, Select Invariant sites for ${1}" >> ${PROGRESSLOG}
}
conda deactivate
conda activate SFStools
# perform the projection; parse the projection and plot the projection
# $1 = input VCF name
# $2 = output file suffix (OUTDIRSNPS previously defined)
easySFS_projection_popHet75() {
# generate easySFS preview
# -a keep all snps
# -v verbose
echo -e "[$(date "+%Y-%m-%d %T")] Start easySFS, perPop maxHet 0.75 for ${1}" >> ${PROGRESSLOG}
python ${WORKSCRIPT} -i ${1} -maxHetFilter ${maxHetFilter} -p ${POPFILE} --preview -a -v > ${OUTDIRSNPS}/PreviewProjection_${2}.txt 2>> ${LOG}

# parse the preview (the population had to match presupplied variables)
INFILE=${OUTDIRSNPS}/PreviewProjection_${2}.txt
echo -e "[$(date "+%Y-%m-%d %T")] Parsing easySFS output ${INFILE}" >> ${PROGRESSLOG}
for pop in ${POPS[@]}
do
OUTFILE=${OUTDIRSNPS}/${pop}_${2}_PreviewProjection_Rformat.txt
echo "projection,snps" > ${OUTFILE}
grep -A1 "$pop$" ${INFILE} | tail -1 | \
sed 's/(//g' | sed 's/)//g' | sed 's/, /,/g' |  tr '\t' '\n' >> ${OUTFILE}
done

echo -e "[$(date "+%Y-%m-%d %T")] Finish easySFS, perPop maxHet 0.75 for ${1}" >> ${PROGRESSLOG}

# perform plotting (redirect stderr to log as well)
#echo -e "[$(date "+%Y-%m-%d %T")] plotting ${INFILE}" >> ${PROGRESSLOG}
#cd ${OUTDIRSNPS}
#Rscript --vanilla ${PLOTSCRIPT} "${OUTDIRSNPS}" "${OUTDIRSNPS}/OptimalProjectionPlots/" ${2} >> ${LOG} 2>&1
}
conda deactivate
###########################################################
## def variables
HOMEDIR=/data/shared/snigenda/finwhale_projects
VCFDIR=${HOMEDIR}/fin_genomics/demographic_inference/Neutral_regions/neutralVCFs
WORKDIR=${HOMEDIR}/fin_genomics/demographic_inference/SFS
WORKSCRIPT=${HOMEDIR}/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/easySFS_a.py
PLOTSCRIPT=${HOMEDIR}/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/easySFS_function_determineOptimalProjection.R
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})

# BLUE WHALE REFERENCE
REFERENCE=${HOMEDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta
EXCLUDE_SAMPLE=("ENPOR12" "ENPCA01" "ENPCA09" "GOC010" "GOC080" "GOC111" "GOC077" "GOC053")
POPS=("ENP" "GOC" "ESP")

OUTDIRSNPS=${WORKDIR}/Neutral/SNPs
OUTDIRINVS=${WORKDIR}/Neutral/Invariants
mkdir -p ${OUTDIRSNPS}
mkdir -p ${OUTDIRINVS}
mkdir -p ${WORKDIR}/Neutral/logs

POPFILE=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/popmap.txt
maxHetFilter=0.75 # still need max het filter if using only PASS sites in filtered vcf, because of population specific het (popHet) could still exceed maxHet 

# additional arguments
mydate=$(date "+%Y%m%d")
LOG=${WORKDIR}/Neutral/logs/easySFS_ProjectionPreview_onlyPASSpopHet75_${mydate}_${IDX}.log
PROGRESSLOG=${WORKDIR}/Neutral/logs/ProjectionPreview_${mydate}_${IDX}_progress.log
###########################################################

## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}; easySFS projection Final (onlyPASS, popHet=0.75)" 

cd ${OUTDIRSNPS}

conda activate whalinGATK
gatk_select_variants_onlyPASS ${VCFDIR}/Neutral_sites_SFS_${IDX}.vcf.gz ${OUTDIRSNPS}/SNPs_neutral_for_SFS_${IDX}.vcf.gz 
conda deactivate

conda activate SFStools
easySFS_projection_popHet75 ${OUTDIRSNPS}/SNPs_neutral_for_SFS_${IDX}.vcf.gz "neutral_SNPS_${IDX}"
conda deactivate

cd ${OUTDIRINVS}
conda activate whalinGATK
gatk_select_invariants_onlyPASS ${VCFDIR}/Neutral_sites_SFS_${IDX}.vcf.gz ${OUTDIRINVS}/Invariants_neutral_for_SFS_${IDX}.vcf.gz

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID} ${IDX} Done"
echo -e "[$(date "+%Y-%m-%d %T")] job for ${SLURM_JOB_ID} ${IDX} Done" >> ${PROGRESSLOG}

conda deactivate
