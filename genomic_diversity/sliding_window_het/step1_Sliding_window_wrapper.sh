#!/bin/bash
#SBATCH --job-name=sliding_window_het
#SBATCH --time=70:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20g
#SBATCH --array=01-21
#SBATCH --partition=CPU
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# author -> Kaden Winspear @ CSUSM -> Eastern pacific fin whale project.
# Date: April 2025

# purpose: The purpose of this wrapper is to perform SlidingWindowHet.py on a directory of passing VCF files. Script calculates the sliding window het within 1mb windows.Tab delimited
#config with scaffold / chromosome names and length is needed. Python 2.7.8 is needed to run the script.

# loading conda environment w/ python 2.7.18 & pysam
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate /home/winsp001/.conda/envs/py27

# Directory to VCFs filtered for PASS / Python script dir / and scaffold lenghts config. To Run config needs to be a two column tab delimineted file consisting of Chr / Scaffold name and length.

VCF_DIR="/data/shared/snigenda/finwhale_projects/filteredvcf/passing_bisnps/passing_vcfs"
SCRIPT_DIR="/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/genomic_diversity/sliding_window_het"
CHROM_LENGTHS="/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/scaffold_lengths.txt"
OUTPUT_DIR="/data/shared/snigenda/finwhale_projects/fin_genomics/diversity/sliding_window/1mbwins"
WINDOW_SIZE=1000000
STEP_SIZE=1000000

# Define VCF files.
VCF_FILE="${VCF_DIR}/JointCalls_08_b_VariantFiltration_$(printf "%02d" $SLURM_ARRAY_TASK_ID)_passing.vcf.gz"

# Check for VCF files.
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF file $VCF_FILE does not exist, uh oh"
    exit 1
fi

# Identify chromosome name from scaffold_len.txt - script uses  NR==ArrayId to find the column that has the chrom / scaffold ID for the array.
CHROM=$(awk "NR==$SLURM_ARRAY_TASK_ID {print \$1}" $CHROM_LENGTHS)

# fails if chrom ID is not found.
if [[ -z "$CHROM" ]]; then
    echo "Error: No chromosome found for SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID , uh oh."
    exit 1
fi

# Runs script.
python $SCRIPT_DIR/step1_SlidingWindowHet.py $VCF_FILE $CHROM_LENGTHS $WINDOW_SIZE $STEP_SIZE $CHROM

