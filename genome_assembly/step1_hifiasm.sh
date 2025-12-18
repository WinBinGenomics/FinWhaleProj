#!/bin/bash
#SBATCH --job-name=hifiasm_ont
#SBATCH --output=hifiasm_ont_%j.out
#SBATCH --error=hifiasm_ont_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=HIMEM
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

#### Load Conda ####
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

#### Load Directories + Reads #####
ont_reads="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/raw_data/finwhale_all.fastq.gz"
outdir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/ont_assembly/"
prefix="gulf_hifiasm_assembly"

mkdir -p "$outdir"

# Runs hifiasm
hifiasm -o "$outdir/$prefix" -t $SLURM_CPUS_PER_TASK --ont "$ont_reads"

echo "Done w/ hifiasm --ont"

conda deactivate
