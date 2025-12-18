#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --output=pixy_%j.out
#SBATCH --error=pixy_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=GPU
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Snigenda workshop
# Author:Kaden Winspear - Dec. 2025
# Usage:sbatch step1_pixy.sh

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate /home/winsp001/.conda/envs/pixy

HOMEDIR="/data/shared/snigenda/finwhale_projects/filteredvcf/workshop"
VCFDIR="${HOMEDIR}/Diversity_analysis/pixy/vcf"
VCF="${VCFDIR}/Pixy_filtered.vcf.gz"
POP_FILE="${HOMEDIR}/scripts/configs/popmap.txt"
OUTPUT_DIR="${HOMEDIR}/Diversity_analysis/pixy"

# Run Pixy [10000 window]
pixy --stats pi --vcf $VCF --populations $POP_FILE --window_size 10000 --n_cores $SLURM_CPUS_PER_TASK --output_folder $OUTPUT_DIR --output_prefix "Chrom19thru21"
