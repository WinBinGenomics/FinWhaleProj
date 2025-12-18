#!/bin/bash
#SBATCH --job-name=currentNe2_array
#SBATCH --output=currentNe2_%x_%A_%a.out
#SBATCH --error=currentNe2_%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --partition=HIMEM
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G


#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Author: Kaden Winspear - Dec. 2025
# Usage: sbatch step5_currentNe2Runner.sh [Population]
# Example: sbatch step5_currentNe2Runner.sh ENP
# Purpose: Runs currentNe2 on populations.

set -euo pipefail

# ------------------------------------------------------
# Add currentNe2 binary to PATH (this job only)
#   (update this path if your CURRENTNE2 dir is named differently)
# ------------------------------------------------------
export PATH=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/programs/currentNe2:$PATH

# ------------------------------------------------------
# Population name (passed as argument)
# ------------------------------------------------------
pop="$1"

# ------------------------------------------------------
# Replicate ID from Slurm array
# ------------------------------------------------------
#run_id=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")

# ------------------------------------------------------
# Input PLINK directory
# ------------------------------------------------------
inputdir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/plink/${pop}"

# ------------------------------------------------------
# Output directory
# ------------------------------------------------------
outputdir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/LD_methods/currentNe2_results/${pop}"
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
outprefix="${outputdir}/${pop}_currentNe_metapop_500k.txt"

# ------------------------------------------------------
# Info logging
# ------------------------------------------------------
echo "Running currentNe2"
echo "Population:      ${pop}"
echo "PED file:        ${pedfile}"
echo "Output prefix:   ${outprefix}"
echo "Threads:         ${SLURM_CPUS_PER_TASK:-1}"
echo "currentNe2 path: $(which currentNe2 || echo 'not found in PATH')"

# ------------------------------------------------------
# Run currentNe2 (mirroring gone2 options)
#   Adjust flags here if currentNe2 uses a different syntax.
# ------------------------------------------------------

currentne2 -t "${SLURM_CPUS_PER_TASK:-1}" "${pedfile}" -r 1.1 -x -s 100000 -o "${outprefix}"
#currentNe2 -t "${SLURM_CPUS_PER_TASK:-1}" "${pedfile}" -r 1.1 -x -s 100000 -S "${seed}" -l 0.002 -u 0.05 -o "${outprefix}"
