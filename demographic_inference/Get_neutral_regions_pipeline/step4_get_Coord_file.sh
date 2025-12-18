#!/bin/bash
#SBATCH --job-name=getNeutralCoords
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --array=1-21
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu


# Gets neutral coordinates to bed file
# Adapted from Annabel Beichman's script (to analyze exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: sbatch get_Coord_file.sh


########## Setting environment
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate py27

set -oe pipefail

########## Set variables, files and directories

HOMEDIR=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference
WORKDIR=${HOMEDIR}/Neutral_regions
VCFDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
OUTDIR=${WORKDIR}/neutralVCFs
SCRIPTDIR=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/Get_neutral_regions_pipeline
NeutralCoord_SCRIPT=${SCRIPTDIR}/step4_obtain_noCpG_noRepetitive_coordinates.py
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})

##### make directories were this information will be stored

mkdir -p ${WORKDIR}/repeatRegions
mkdir -p ${WORKDIR}/get_fasta


##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

echo "[$(date "+%Y-%m-%d %T")] Start creating bed files for ${SLURM_ARRAY_TASK_ID} JOB_ID: ${SLURM_JOB_ID}"
echo "The sbatch input"
echo "${SLURM_ARRAY_TASK_ID}"

PROGRESSLOG=./logs/create_beds_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${SLURM_JOB_ID}" > ${PROGRESSLOG}


########## Creates a bed file for each vcf file

echo -e "[$(date "+%Y-%m-%d %T")]  creating bed files ... " >> ${PROGRESSLOG}
LOG=./logs/Step2_Creating_bed_files_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

python ${NeutralCoord_SCRIPT} --VCF ${VCFDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz --outfile ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_${IDX}.bed

conda deactivate

########## deactivate py2.7 conda env & Activates conda w/ bedtools
conda activate fintools

bedtools merge -d 10 -i ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_${IDX}.bed > ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_mergedMaxDistance10_${IDX}.bed

cat ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_mergedMaxDistance10_*.bed | sort -k1,1 -k2,2n > ${WORKDIR}/get_fasta/min20kb_DistFromCDRs_noCpGIsland_noRepeat_mergedMaxDistance10_sorted.bed 


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


conda deactivate

