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
conda activate fintools

set -o pipefail


#Defining directories ---------------------

workdir=/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH


#Main ---------------------------------------

for pop in GOC ENP ESP
do
	#New directory
	outdir=$workdir/$pop
  	cd ${outdir}

	#Get all the vcf files made in the previous step
	ls *.vcf.gz > vcf_files.txt

	#Concat vcf files
 	bcftools concat -f vcf_files.txt -O z -o JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz
	tabix -p vcf JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz

	#Convert concat vcf file into Oxford file
	bcftools convert \
	--no-version \
	-g JointCalls_08_B_VariantPASS_${pop}_allcontig \
	--3N6 \
	--vcf-ids \
	--threads 8 \
	JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz
done

conda deactivate
