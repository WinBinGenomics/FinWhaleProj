#!/bin/bash
#SBATCH --job-name=easySFSProjection_monomorphic
#SBATCH --time=23:00:00
#SBATCH --partition=HIMEM
#SBATCH --mem=15G
#SBATCH --array=1-21
#SBATCH --output=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/reports/SFS_monomorphic_%A_%>
#SBATCH --error=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/reports/SFS_monomorphic_%A_%a>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --cpus-per-task=1

# This script calculates monomorphic site.
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Adjsuted by kaden Winspear @ CSUSM for SLURM server. Reminder -> script is written in 2.7.8
# Usage: sbatch SFS_add_monomorphic_sites.sh

module load snigenda-anaconda3/2024.02
source /cm/shared/apps/snigenda-anaconda3/2024.02/etc/profile.d/conda.sh
conda activate py27


set -o pipefail

#Define directories ---------------------------------

homedir=/data/shared/snigenda/finwhale_projects
workdir=${homedir}/fin_genomics/demographic_inference/SFS/Neutral
vcfdir=${homedir}/fin_genomics/demographic_inference/Neutral_regions/neutralVCFs
outdir=${workdir}/SFS_projection_Monomorphic
#monoscript=${homedir}/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/getMonomorphicProjectionCounts.1D.2DSFS.py
monoscript=${homedir}/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/getMonomorphicProjectionCounts.1D.2D.3D.py # modification to count 3D monosites
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})


#Main ---------------

popfile=${homedir}/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/popmap.txt
projections="46,34,26"

mkdir -p ${outdir}
cd ${outdir}

########## get counts of monomorphic sites to add to the SFSes ############

python ${monoscript} --vcf ${vcfdir}/Neutral_sites_SFS_${IDX}.vcf.gz --popMap ${popfile} --proj ${projections} --popIDs ENP,ESP,GOC --outdir ${outdir} --outPREFIX ${IDX}

