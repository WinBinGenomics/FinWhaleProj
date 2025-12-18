#!/bin/bash
#SBATCH --account=winsp001
#SBATCH --job-name=pixy_array
#SBATCH --ntasks-per-node=1
#SBATCH --array=01-21
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=22:00:00
#SBATCH --partition=CPU
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# author: Kaden Winspear - > CSUSM. Apr. 04 2025
# USAGE: sbatch pixy.sh | uses pixy to estimate genomewide pi | dxy | fst values for the Eastern Pacific Fin Whale Populations. Window size is largest scaffold / chrom.


module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate /home/winsp001/.conda/envs/pixy

VCFDIR="/data/users/snigenda/finwhale_data/blue_ref_omar_data/blue_ref_vcfs"

# Get VCFs to loop over.
IDX=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
VCF="${VCFDIR}/JointCalls_08_b_VariantFiltration_${IDX}.vcf.gz"

# Popmap -> pixy input
POP_FILE="/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/pixy_popmap.txt"

# output directory
OUTPUT_DIR="/data/shared/snigenda/finwhale_projects/fin_genomics/diversity/genomewide_diversity/pixy_estimates/genomewide"

# Runs Pixy
pixy --stats pi fst dxy \
     --vcf $VCF \
     --populations $POP_FILE \
     --window_size 18500000 \
     --n_cores $SLURM_CPUS_PER_TASK \
     --output_folder $OUTPUT_DIR \
     --output_prefix "Whole_chrom_Analysis_${SLURM_ARRAY_TASK_ID}" \
