#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH --output=medaka_%j.out
#SBATCH --error=medaka_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=HIMEM
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

#### Load Conda ####
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

# Path to respective DIRs
assemblyDir='/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/ont_assembly'
ont_reads="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/raw_data/finwhale_all.fastq.gz"
assembly="${assemblyDir}/gulf_hifiasm_assembly.fa"
outdir="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/medaka"

mkdir -p "$outdir"

# Convert GFA to FASTA
awk '/^S/{print ">"$2"\n"$3}' "${assemblyDir}/gulf_hifiasm_assembly.bp.p_ctg.gfa" > "$assembly"

# Run Medaka
medaka_consensus \
  -i "$ont_reads" \
  -d "$assembly" \
  -o "$outdir" \
  -t $SLURM_CPUS_PER_TASK

if [[ $? -eq 0 ]]; then
    echo "IT WORKED"
else
    echo "Medaka failed" >&2
    exit 1
fi

conda deactivate
