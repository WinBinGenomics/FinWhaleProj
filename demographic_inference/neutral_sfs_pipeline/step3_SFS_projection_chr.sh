#!/bin/bash
#SBATCH --job-name=easySFSProjection_chr
#SBATCH --time=05:00:00
#SBATCH --mem=5G
#SBATCH --output=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/reports/SFS_projection_%A_%a.out.txt
#SBATCH --error=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/reports/SFS_projection_%A_%a.err.txt
#SBATCH --array=1-21
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu


# This script runs SFS projection  per chromosomes
# Author: Meixi Lin, modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Adjusted by Kaden Winspear @ CSUSM for SLURM. Eastern pacific fin whale project.
# Usage: sbatch SFS_projection_chr.sh

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate SFStools

set -o pipefail

#Define directories ---------------------------------
homedir=/data/shared/snigenda/finwhale_projects
workdir=${homedir}/fin_genomics/demographic_inference/SFS/Neutral
scriptdir=${homedir}/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline
easySFS=${scriptdir}/easySFS_a.py
vcfdir=${workdir}/SNPs
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})


#Main ---------------

cd ${workdir}
#popfile=${scriptdir}/popmap.txt
popfile=${scriptdir}/pop_map.txt
outdir=${workdir}/SFS_projection

mkdir -p $outdir
cd ${outdir}

#python $easySFS -i ${vcfdir}/SNPs_neutral_for_SFS_${IDX}.vcf.gz -p ${popfile} --proj 46,34,26 -a -f -v -o SNPs_easySFS_projection_${IDX} -maxHetFilter 0.75
python $easySFS -i ${vcfdir}/SNPs_neutral_for_SFS_${IDX}.vcf.gz -p ${popfile} --proj 34,46,26 -a -f -v -o SNPs_easySFS_projection_${IDX} -maxHetFilter 0.75

conda deactivate

