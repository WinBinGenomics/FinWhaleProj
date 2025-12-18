#!/bin/bash
#SBATCH --job-name=gone2_array
#SBATCH --output=gone2_%x_%A_%a.out
#SBATCH --error=gone2_%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --partition=HIMEM
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-100


#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Author: Kaden Winspear - Dec. 2025
# Usage: sbatch step3_goneRunner.sh [Population] Example: sbatch step3_goneRunner.sh ENP
# Purpose: Runs gone2 on populations.

set -euo pipefail

# ------------------------------------------------------
# Add gone2 binary to PATH (this job only)
# ------------------------------------------------------
export PATH=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/programs/GONE2:$PATH

# ------------------------------------------------------
# Population name (passed as argument)
# ------------------------------------------------------
pop="$1"

# ------------------------------------------------------
# Replicate ID from Slurm array
# ------------------------------------------------------
run_id=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")

# ------------------------------------------------------
# Input PLINK directory
# ------------------------------------------------------
inputdir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/plink/${pop}"

# ------------------------------------------------------
# Output directory
# ------------------------------------------------------
outputdir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/gone_results/${pop}"
mkdir -p "${outputdir}"

# ------------------------------------------------------
# Input PED file
# ------------------------------------------------------
pedfile="${inputdir}/${pop}_autosomes_LD.ped"

if [[ ! -f "${pedfile}" ]]; then
    echo "ERROR: PED file not found: ${pedfile}" >&2
    exit 1
fi

# ------------------------------------------------------
# Per-run output prefix
# ------------------------------------------------------
outprefix="${outputdir}/${pop}_run${run_id}"

# ------------------------------------------------------
# Info logging
# ------------------------------------------------------
echo "Running GONE2"
echo "Population:   ${pop}"
echo "Run ID:       ${run_id}"
echo "PED file:     ${pedfile}"
echo "Output:       ${outprefix}"
echo "Threads:      ${SLURM_CPUS_PER_TASK:-1}"
echo "gone2 path:   $(which gone2)"

# ------------------------------------------------------
# Run GONE2 (with randomness via -x)
# ------------------------------------------------------
seed=$((10000 + SLURM_ARRAY_TASK_ID))

gone2 -t "${SLURM_CPUS_PER_TASK:-1}" "${pedfile}" -r 1.1 -x  -s 100000 -S "${seed}" -o "${outprefix}"
#gone2 -t "${SLURM_CPUS_PER_TASK:-1}" "${pedfile}" -r 1.1  -s 100000 -S "${seed}" -u 0.03 -l 0.001 -o "${outprefix}"
