#!/bin/bash
#SBATCH --job-name=fastqc_check
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=HIMEM
#SBATCH --mem=400G
#SBATCH --mail-user=winsp001@csusm.com
#SBATCH --mail-type=ALL

set -euo pipefail

#### Load Conda ####
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

# Input directory and files
#indir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/illumina_data/trimmed"
indir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/raw_data"
#raw_R1="SH5556_R1_paired.fastq.gz"
#raw_R2="SH5556_R2_paired.fastq.gz"
ontdata="finwhale_all.fastq.gz"

# Create an output folder
outdir="${indir}/fastqc_reports"
mkdir -p "$outdir"

# Run FastQC on both files
#fastqc "${indir}/${raw_R1}" "${indir}/${raw_R2}" --outdir "$outdir"
fastqc "${indir}/${ontdata}" --nano --outdir "$outdir"


conda deactivate
