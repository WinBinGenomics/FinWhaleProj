#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --output=quast_%j.out
#SBATCH --error=quast_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

#### Load Conda ####
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

# --- Paths ---
OUTDIR="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/quast"
mkdir -p "$OUTDIR"

ASM1="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/ont_assembly/gulf_hifiasm_assembly.fa"
ASM2="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/medaka/consensus.fasta"
ASM3="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/np2/gulf.np2.fa"
ASM4="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/ont_assembly/40-polishing/filtered_contigs.fasta"
# --- Run QUAST ---
quast.py -t $SLURM_CPUS_PER_TASK \
  -o "$OUTDIR" \
  "$ASM1" "$ASM2" "$ASM3" "$ASM4"

echo "QUAST finished. Results in $OUTDIR"
