#!/bin/bash
#SBATCH --mem=40G
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

module load raxml/8.2.13

phyldir='/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/all70/Blue/raxml_data'
phyfile="${phyldir}/JointCalls_all70_filterpass_bialleic_all_LDPruned_maf05_SA_mrF.min4.phy"

# Run RAxML with 100 bootstraps
raxmlHPC -s ${phyfile} -n snps_maf05_LDpruned -m GTRGAMMA -p 12345 -T 8 -N 100 -w ${phyldir}
