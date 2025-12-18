#!/bin/bash
#SBATCH --job-name=fscWrapper3D_SplitNomig
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --partition=GPU
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-100

# Wrapper to run 100 fastsimcoal iterations for a 3D joint SFS model

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate SFStools

set -o pipefail

# Defined directories -------------
wd=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal
infDir=$wd/                                   # where inference is happening
genericDir=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/fsc26/3DModels
sfs=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection_fsc/ESP-ENP-GOC.MSFS.obs

# Parameters -----------------------
models='Split-symm-m01-asym-m12'
cores=8

########################################  MAIN  #############################################

for model in $models
do
    echo "starting $model"
    header=${model}

    # Copy generic files into directory and update
    outdir=$infDir/3DModels/$model/run_${SLURM_ARRAY_TASK_ID}
    mkdir -p "$outdir"

    cp "$genericDir/$model.tpl" "$genericDir/$model.est" "$outdir"   # copy .est and .tpl files to outdir

    # Copy SFS into inference directory and rename to match .tpl/.est header
    cp "$sfs" "$outdir/${header}_MSFS.obs"

    cd "$outdir"

    # Run fastsimcoal2
    #fsc27093  -t "${header}.tpl" -e "${header}.est" -n 1000000 -m -M -L 60 -c${cores} -q
    fsc27093 -t ${header}.tpl -n 1000000 -m --multiSFS -e ${header}.est -M -L 60 -c${cores} -q
done


