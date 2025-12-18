#!/bin/bash
#SBATCH --job-name=yak_count
#SBATCH --output=yak_count_%j.out
#SBATCH --error=yak_count_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=HIMEM
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G

set -euo pipefail

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

# Paths
trimdir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/illumina_data/trimmed"
outdir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/yak"

mkdir -p "$outdir"

# Input paired reads
R1="$trimdir/SH5556_R1_paired.fastq.gz"
R2="$trimdir/SH5556_R2_paired.fastq.gz"

# Output
yakdb="$outdir/illumina.yak"

# Run yak count
yak count -k31 -b37 -t $SLURM_CPUS_PER_TASK -o "$yakdb" \
    <(zcat "$R1") \
    <(zcat "$R2")

echo "Yak count finished. Output: $yakdb"
