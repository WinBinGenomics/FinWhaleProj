#! /bin/bash
#SBATCH --mem=24G
#SBATCH --time=23:00:00
#SBATCH --partition=HIMEM
#SBATCH --array=01-21
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=winsp001@csusm.edu


# Author: Sergio Nigenda, based on Tanya Phung's scripts for step 10 of her NGS pipeline.
# Adjusted / adapted by Kaden Winspear @ CSUSM - > Eastern Pacific Fin Whale project.

########## Purpose
# Generates a bed file with the high quality sites in the VCF files
# Usage example: qsub GetHiQualCoords_20200602.sh Blue


########## Setting conda environment

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate py27

set -o pipefail

########## Set variables, directories and files

#REF=${1}

HOMEDIR=/data/shared/snigenda/finwhale_projects
WORKDIR=${HOMEDIR}/fin_genomics/demographic_inference/Neutral_regions/HiQualCoords
SCRIPTDIR=${HOMEDIR}/scripts/winfingenomics/demographic_inference/Get_neutral_regions_pipeline
VCFDIR=/data/users/snigenda/finwhale_data/blue_ref_omar_data/blue_ref_vcfs
HiQualCoords_SCRIPT=${SCRIPTDIR}/obtain_high_qual_coordinates.py
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})

mkdir -p ${WORKDIR}

########## Logging

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Start GetHiQCoords for ${REF} {IDX} Job ID: ${SLURM_ARRAY_TASK_ID}"
echo "The sbatch input"
echo "${SLURM_ARRAY_TASK_ID}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/GetHiQualSitesCoords_A_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${SLURM_ARRAY_TASK_ID}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] Selecting high quality coordinates ... " >> ${PROGRESSLOG}

LOG=./logs/01_A_Get_HighQuality_Coords_${REF}_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}


########## Step 9A. Selecting high quality coordinates and missing sites. If ${HiQualCoords_SCRIPT} is defined as obtain_high_qual_coordinates.py it selects only high quality data to later build neutral SFS without projection. However, if is defined as obtain_high_qual_coordinates_miss.py then it includes high quality and missing data for SFS projection. 

python ${HiQualCoords_SCRIPT} --VCF ${VCFDIR}/JointCalls_08_b_VariantFiltration_${IDX}.vcf.gz --outfile HQsitesCoords_${IDX}.bed

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

conda deactivate

########## Step 9B. Merging High quality coordinates (The previous step is made for each variant, therefore here we merge those coordinates to have less individual regions)

conda activate fintools

PROGRESSLOG=./logs/GetHiQualSitesCoords_B_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] Merging high quality coordinates ... " >> ${PROGRESSLOG}

LOG=./logs/01_B_Merge_HighQuality_Coords_${REF}_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

bedtools merge -i HQsitesCoords_${IDX}.bed > HQsitesCoords_merged_${IDX}.bed

wait

cat HQsitesCoords_merged_*.bed | sort -k1,1 -k2,2n > all_HQCoords_sorted_merged.bed

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

conda deactivate
