#!/bin/bash
#SBATCH --job-name=yak_qv
#SBATCH --output=yak_qv_%j.out
#SBATCH --error=yak_qv_%j.err
#SBATCH --partition=HIMEM
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G

set -euo pipefail

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass


# Inputs
yakdb="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/yak/illumina.yak"
asm_hifiasm="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/ont_assembly/gulf_hifiasm_assembly.fa"   # from GFA→FA
asm_medaka="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/medaka/consensus.fasta"

# Outputs
outdir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/medaka/qv"
mkdir -p "$outdir"

# QV before Medaka (hifiasm draft)
yak qv -t $SLURM_CPUS_PER_TASK -p -K3g -l100k "$yakdb" "$asm_hifiasm" > "$outdir/qv.hifiasm.txt"

# QV after Medaka
yak qv -t $SLURM_CPUS_PER_TASK -p -K3g -l100k "$yakdb" "$asm_medaka"  > "$outdir/qv.medaka.txt"

if [[ $? -eq 0 ]]; then
    echo "IT WORKED"
else
    echo "failed :(" >&2
    exit 1
fi

conda deactivate
