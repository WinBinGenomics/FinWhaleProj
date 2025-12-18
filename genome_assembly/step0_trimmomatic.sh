#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --output=trim_%j.out
#SBATCH --error=trim_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=GPU
#SBATCH --mem=16G
#SBATCH --mail-user=winsp001@csusm.com
#SBATCH --mail-type=ALL

#### Load Conda ####
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

#### Load respective directories ####
datadir='/data/users/snigenda/finwhale_data/finwhale_raw_gnom_seq/FT-SH5556/FT-SA39975-FT-SPN00435_HTJWGDSXX'
outdir='/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/illumina_data/trimmed'
adapters='/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/config/TruSeq-PE.fa'

raw_R1='SH5556_SA39975_S9_L003_R1_001.fastq.gz'
raw_R2='SH5556_SA39975_S9_L003_R2_001.fastq.gz'


# Ensure outdir exists
mkdir -p "$outdir"

# Go to data dir
cd "$datadir"

# Run Trimmomatic
trimmomatic PE -threads 8 \
  "$raw_R1" "$raw_R2" \
  "$outdir/SH5556_R1_paired.fastq.gz" "$outdir/SH5556_R1_unpaired.fastq.gz" \
  "$outdir/SH5556_R2_paired.fastq.gz" "$outdir/SH5556_R2_unpaired.fastq.gz" \
  ILLUMINACLIP:"$adapters":2:30:10:8:true \
  SLIDINGWINDOW:4:15 MINLEN:36

conda deactivate
