#!/bin/bash
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --job-name=countsites_py
#SBATCH --mail-type=ALL
#SBATCH --array=1-23%5
#SBATCH --mail-user=winsp001@csusm.edu


# @description	Wrapper of calling the countSitesPerIndividual.py in CSUSM HPC
# Initial script written by the author: Meixi Lin
# Adjusted / Modified by Kaden Winspear @ CSUSM. Apart of Eastern Pacific Fin Whale Project.

#Usage: this is a wrapper for a python script written in py 2.7. To use on slurm cluster -> sbatch step1_wrapper_countSitesPerIndividual.sh

###########################################################
## import packages
#sleep $((RANDOM % 120

module load snigenda-anaconda3/2024.02
source /cm/shared/apps/snigenda-anaconda3/2024.02/etc/profile.d/conda.sh
conda activate /cm/shared/apps/snigenda-anaconda3/2024.02/envs/py2.7.18

set -eo pipefail

###########################################################
## def functions

###########################################################
## input variables
DATASET="all70"
REF="Blue"

## def variables
TODAY=$(date "+%Y%m%d")
WORK="/data/shared/snigenda/finwhale_projects/fin_genomics/heterozygosity/"
WORKSCRIPT="/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/genomic_diversity/genome_stats/genomewide_diversity/step1_countSitesPerIndividual.py"
OUTDIR="${WORK}/genomewide_diversity/count_sites_${TODAY}"
VCFDIR="/data/users/snigenda/finwhale_data/blue_ref_omar_data/blue_ref_vcfs"

IDX=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
VCF="JointCalls_08_b_VariantFiltration_${IDX}.vcf.gz"

mkdir -p ${OUTDIR}

###########################################################
## main
# first check the overall vcf.gz file
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} Count sites for ${VCF}"

python ${WORKSCRIPT} --vcf ${VCFDIR}/${VCF} --outfile ${OUTDIR}/${DATASET}_${REF}_${IDX}_sites_summary.txt --filter "PASS" --contig ${IDX}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SLURM_ARRAY_TASK_ID} Done"

conda deactivate
