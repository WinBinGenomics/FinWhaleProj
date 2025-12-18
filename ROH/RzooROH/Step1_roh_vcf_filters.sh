#! /bin/bash
#SBATCH --job-name=filter_vcftools
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=HIMEM
#SBATCH --array=1-21
#SBATCH --mail-user=winsp001@csusm.edu


# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Adjusted / Modified by Kaden Winspear @ CSUSM. Eastern Pacific fin whale project.
# Date: May 1st 2025
# Purpose: This script filters vcf and subset them per populations.

module load vcftools/0.1.16

set -o pipefail

#Defining directories ---------------------

workdir=/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH
vcfdir=/data/shared/snigenda/finwhale_projects/filteredvcf/passing_bisnps
IDX=$(printf %02d ${SLURM_ARRAY_TASK_ID})


#Main ---------------------------------------

for pop in GOC ENP ESP
do
  outdir=$workdir/$pop
  mkdir -p $outdir
  cd ${outdir}

#adjusted vcftools command for biallelic snps only.

vcftools --gzvcf ${vcfdir}/JointCalls_08_b_VariantFiltration_${IDX}_passing_bisnps.vcf.gz \
  --remove-filtered-all \
  --keep ${pop}_samples.txt \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out ${workdir}/${pop}/VariantFiltration_${IDX}_OnlyPass.${pop}

  gzip ${workdir}/${pop}/VariantFiltration_${IDX}_OnlyPass.${pop}.recode.vcf

done

