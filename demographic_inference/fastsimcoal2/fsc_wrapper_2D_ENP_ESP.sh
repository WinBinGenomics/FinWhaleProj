#!/bin/bash
#SBATCH --job-name=fscWrapper2D
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --partition=GPU
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-100


# This is a wrapper that will run 100 fastsimcoal iterations for joint SFS for any list of models
# Author: Anabell Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: sbatch fsc_wrapper.sh


module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate SFStools

set -o pipefail


# Defined directories -------------

wd=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal
infDir=$wd/ # specify where inference is happening
genericDir=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/fsc26/2DModels/ENP-ESP
sfs=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection_fsc/ENP-ESP_jointMAFpop1_0.obs

# Parameters -----------------------
pops='ENP-ESP'
models='AncestralSizeChange-Split-iso-assym'
cores=8
version="Normal"

########################################  MAIN  #############################################


for model in $models
do

        echo "starting $model"
        header=${model}

        # Copy generic files into directory and update
        outdir=$infDir/2DModels/$model/$version/run_${SLURM_ARRAY_TASK_ID}
        mkdir -p $outdir

        cp $genericDir/$model.tpl $genericDir/$model.est $outdir # copy .est and .tpl files to outdir

        # Get sfs into inference directory and rename to match .est and .tpl files
        cp $sfs $outdir/${header}_jointMAFpop1_0.obs # copy your sfs into the directory where you'll be doing the fsc inference
        cd $outdir
	fsc27093 -t ${header}.tpl -n 1000000 -m -e ${header}.est -M -L 60 -c${cores} -q

done

