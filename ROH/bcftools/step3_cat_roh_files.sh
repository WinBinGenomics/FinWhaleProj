#!/bin/bash
#SBATCH --job-name=concat_ROH_bcftools
#SBATCH --time=24:00:00
#SBATCH --partition=HIMEM
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu

# This script concatenates all output files per chromosome/scaffold from bcftools roh into a single file
# Author: Sergio Nigenda & Adjusted by Kaden Winspear @ CSUSM -> Eastern Pacific fin whale project.
# Usage: qsub cat_roh_files.sh populations

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

set -o pipefail


########## Defining variables, directories and files

pop=$1

homedir=/data/shared/snigenda/finwhale_projects/fin_genomics
workdir=${homedir}/ROH/rohbcftools
scriptdir=${homedir}/scripts/winfingenomics/ROH/bcftools
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})
popdir=${workdir}/${pop}

cd ${popdir}

cat ${pop}_roh_bcftools_AFACANGT_{01..21} > ${pop}_concat_fwhale_roh_bcftools_G30_ACANGT_1.out

wait

grep -v "^#" ${pop}_concat_fwhale_roh_bcftools_G30_ACANGT_1.out > ${pop}_concat_fwhale_roh_bcftools_G30_ACANGT.out
