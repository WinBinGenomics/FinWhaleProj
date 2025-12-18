#!/bin/bash
#SBATCH --job-name=flye_ont
#SBATCH --output=flye_ont_%j.out
#SBATCH --error=flye_ont_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G

# (use your env that has flye)
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

# Paths
ont_reads="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/raw_data/finwhale_all.fastq.gz"
outdir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/ont_assembly/"
prefix="gulf_flye_assembly"
GENOME_SIZE="2.7g"

mkdir -p "$outdir"

# Flye assembly (use --nano-hq for Dorado sup/duplex reads)
flye \
  --nano-hq "$ont_reads" \
  --genome-size "$GENOME_SIZE" \
  --threads 48 \
  --out-dir "$outdir"

if [[ $? -eq 0 ]]; then
    echo "IT WORKED"
else
    echo "failed :(" >&2
    exit 1
fi

conda deactivate
