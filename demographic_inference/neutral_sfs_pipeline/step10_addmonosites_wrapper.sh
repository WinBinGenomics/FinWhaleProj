#!/bin/bash
#SBATCH --partition=GPU

# wrapper to run the conca.py to add 2Dsfs together.
# usage: sbatch step9_wrapper.sh ENP-GOC or ENP-ESP or ESP-GOC or GOC-ESP

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate SFStools

homedir="/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection/all_dadi_SFS"
script="/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/neutral_sfs_pipeline/add_monosites_3D.py"

inputdir="${homedir}/$1"

if [[ -d "$inputdir" ]]; then
  echo "=== Processing $1 ==="
  cd "$inputdir" || exit 1
  python "$script"
else
  echo "Directory not found: $inputdir"
  exit 1
fi
