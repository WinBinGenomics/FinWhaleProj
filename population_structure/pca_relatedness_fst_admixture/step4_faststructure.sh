#!/bin/bash
#SBATCH --mem=20G
#SBATCH --time=23:00:00
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

## Set environment and handle errors
set -eo pipefail

# Load modules
module load R/4.3.1+Bioconductor
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fastSTRUCTURE

# Basic configuration
DATASET='all70'
REF='Blue'
MAF='05'  # Keep as '05' if your files use this format
WORKDIR="/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/${DATASET}/${REF}"
OUTDIR="${WORKDIR}/faststructure/maf${MAF}"

# Create output directory
mkdir -p "${OUTDIR}"
cd "${OUTDIR}" || exit 1

# Input file check
BED_FILE="${WORKDIR}/JointCalls_${DATASET}_biallelic_all_LDPruned_maf${MAF}_SA_mrf.bed"
if [[ ! -f "$BED_FILE" ]]; then
    echo "ERROR: Missing BED file: $BED_FILE"
    exit 1
fi

# Simple loop structure
for K in {2..6}; do
    for ITER in {1..10}; do
        structure.py \
            -K $K \
            --input="${BED_FILE%.bed}" \
            --cv=5 \
            --seed=$RANDOM \
            --prior=simple \
            --output="K${K}_iter${ITER}"
    done
done

echo "Basic analysis completed"
