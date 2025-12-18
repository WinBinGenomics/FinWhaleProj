#!/bin/bash
#SBATCH --job-name=NoCpGRepeats
#SBATCH --time=23:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --array=03-21

# Extract genomic regions that are at least 10 Kb apart from exons and are not within repetitive regions or cpg islands, to new vcf files per chromosome or scaffold
# Adapted from Annabel Beichman's script (for analyzing exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: sbatch Extract_noCpG_noRepetitive.sh [reference species]

# Edited / modified by Kaden Winspear @ CSUSM for Eastern pacific fin whale project. SLURM cluster.

########## Setting environement ##############

module load snigenda-anaconda3/2024.02
source /cm/shared/apps/snigenda-anaconda3/2024.02/etc/profile.d/conda.sh
conda activate GentoolsPop


set -oe pipefail


########## Set variables, files and directories

REF=${1}

HOMEDIR=/data/shared/snigenda/finwhale_projects
WORKDIR=${HOMEDIR}/fin_genomics/demographic_inference/Neutral_regions
#VCFDIR=${HOMEDIR}/filteredvcf/all70/${REF}
VCFDIR=/data/users/snigenda/finwhale_data/blue_ref_omar_data/blue_ref_vcfs/
OUTDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})
twentyKb=${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromCDR.0based.bed

mkdir -p ${WORKDIR}
mkdir -p ${OUTDIR}


if [ $REF == 'Minke' ]; then
    REFDIR=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0
    REFERENCE=${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic.fasta
    CPG_REPEATS_ALL=${REFDIR}/repeats_cpg_all.bed
fi

if [ $REF == 'Bryde' ]; then
    # Note that for the Bryde's whale, we only used the WMdust.bed file
    REFDIR=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data
    REFERENCE=${REFDIR}/Balaenoptera_edeni_HiC.fasta
    GFF=${REFDIR}/Balaenoptera_edeni/genome/maker/Balaenoptera_edeni_Balaenoptera_edeni_HiC.fasta_v2.functional.gff3.gz
    MASK=${REFDIR}/Balaenoptera_edeni_HiC_repeats/Balaenoptera_edeni_HiC_repeats_WMdust.bed

fi

if [ $REF == 'Blue' ]; then
    REFDIR=${HOMEDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3
    REFERENCE=${REFDIR}/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta
    CPG_REPEATS_ALL=${REFDIR}/MaskFiles/repeats_cpg_all.bed
fi

##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Start extracting nor repetitive regions neither cpg islands for ${REF} ${SLURM_ARRAY_TASK_ID} JOB_ID: ${JOB_ID}"
echo "The SBATCH input"
echo "${REF} ${SLURM_ARRAY_TASK_ID}"

PROGRESSLOG=./logs/Extract_No_repeats_cpg_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${SLURM_JOB_ID}" > ${PROGRESSLOG}


########## Obtain vcf files per chromosome or scaffold that do not contain repeat regions or cpg islands and are at least 10 Kb apart from exons

echo -e "[$(date "+%Y-%m-%d %T")]  Extracting neutral regions with GATK SelecVariants... " >> ${PROGRESSLOG}
LOG=./logs/Step1_ExtractNeutralSites_${REF}_SelectVariants_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk --java-options "-Xmx15g -Djava.io.tmpdir=./temp" SelectVariants \
  -R ${REFERENCE} \
  -V ${VCFDIR}/JointCalls_08_b_VariantFiltration_${IDX}.vcf.gz \
  -XL ${CPG_REPEATS_ALL} \
  -L ${twentyKb} \
  -O ${OUTDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz &>> ${LOG}


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo "[$(date "+%Y-%m-%d %T")] Done extracting no repeats regions or cpg islands for ${REF} ${SLURM_ARRAY_TASK_ID} Job ID: ${SLURM_JOB_ID}"

conda deactivate
