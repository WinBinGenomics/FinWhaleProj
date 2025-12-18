#!/bin/bash
#SBATCH --partition=CPU

# This script is used add all the 2D and 3D sfs per chromosomes. And places them into a specific directory (Either ENP-GOC, ENP-ESP, or ESP-GOC).

set -euo pipefail

base="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection"
outdir="${base}/all_dadi_SFS"
mkdir -p "$outdir"

for i in $(seq -w 01 21); do
  dir="${base}/SNPs_easySFS_projection_${i}/dadi"
  sfs="${dir}/ESP-ENP-GOC.sfs" # adjust as needed.

  if [ -f "$sfs" ]; then
    target="${outdir}/ESP-ENP-GOC_chr${i}.sfs" # adjust as needed
    echo "Copying: $sfs -> $target"
    cp -- "$sfs" "$target"
  else
    echo "No sfs in $dir"
  fi
done
