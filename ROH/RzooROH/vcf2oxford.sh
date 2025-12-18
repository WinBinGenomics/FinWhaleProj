#! /bin/bash
#SBATCH --job-name=concat2oxford
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --partition=HIMEM
#SBATCH --mail-user=winsp001@csusm.edu

# Author: Meixi Lin, modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Adjusted by Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale Project
# Date: May 1st 2025

# Purpose: This script concat vcf and convert to oxford files


module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate bcftools_v1.9


set -o pipefail


#Defining directories ---------------------

workdir=/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH


#Main ---------------------------------------

for pop in GOC ENP ESP
do
    outdir=$workdir/$pop
    cd "${outdir}"

	bcftools convert \
	-g JointCalls_08_B_VariantPASS_${pop}_allcontig \
	--chrom --tag GT \
	--threads 8 \
	JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz

done

conda deactivate
