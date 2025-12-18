#!/bin/bash
#SBATCH --mem=40G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

vcfdir='/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/all70/Blue/'
outputdir='/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/all70/Blue/raxml_data/'
vcffile="${vcfdir}/JointCalls_all70_filterpass_bialleic_all_LDPruned_maf05_SA_mrF.vcf"
scriptdir='/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/population_structure/rax_ml_tree/vcf2phylip/'

# Ensure output directory exists
mkdir -p "${outputdir}"

"${scriptdir}/vcf2phylip.py" -i "${vcffile}" -o "snps_maf05_LDpruned"
