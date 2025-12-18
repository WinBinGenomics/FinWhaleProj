#!/bin/bash
#SBATCH --job-name=easySFSProjection
#SBATCH --time=05:00:00
#SBATCH --mem=15G
#SBATCH --output=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/reports/SFS_projection_%A_%a.out.txt
#SBATCH --error=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/reports/SFS_projection_%A_%a.err.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu


# This script concentrate files names of SFS projection per chromosomes
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Adjusted by Kaden Winspear @ CSUSM for SLURM. Eastern pacific fin whale project.
# Usage: sbatch SFS_projection_join.sh

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate SFStools

set -o pipefail

#Define directories ---------------------------------
workdir=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection_fsc
SFSdir=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection

mkdir -p ${workdir}
mkdir -p ${SFSdir}

#Main ---------------

for pop in GOC ENP ESP
do
	cd ${workdir}

  	x=$(ls $SFSdir) ; printf "%s\n" "$x" > directories.txt
	while read line; do
		echo -e $SFSdir"/"$line"/fastsimcoal2/"$pop"_MAFpop0.obs" >> $pop"_SFS_projection_files.txt"
	done < directories.txt
done

conda deactivate
