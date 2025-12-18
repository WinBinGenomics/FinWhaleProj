#!/bin/bash
#SBATCH --job-name=fscWrapperGOC
#SBATCH --output=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/fsc26/1DModels/reports/fscwrapperGOC.%A_%a.out.txt
#SBATCH --error=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/fsc26/1DModels/reports/fscwrapperGOC.%A_%a.err.txt
#SBATCH --time=10:00:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=8
#SBATCH --partition=GPU
#SBATCH --mail-user=winsp001@hpc-login.csusm.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-100


# This is a wrapper that will run 100 fastsimcoal iterations for each population for any list of models
# Author: Anabell Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: sbatch fsc_wrapper_GOC.sh
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate SFStools


set -o pipefail


# Defined directories -------------

wd=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal
infDir=$wd/1DModels # specify where inference is happening
genericDir=/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/demographic_inference/fsc26/1DModels # location of generic FSC models
sfsDir=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection_fsc #sfsDir=/u/home/p/pnunez/project-rwayne/FinWhale/fastsimcoal/optimization_tryout

# Programs -------------------------
#fsc=/u/home/p/pnunez/programs/fsc26_linux64/fsc26

# Parameters -----------------------

models='1D.4Epoch.double'
pops="GOC"
cores=8
ss="26"
version="Normal"

########################################  MAIN  #############################################


for pop in $pops
do

	# Get sample size for the population

	for model in $models
	do

		echo "starting $pop, $model"
		header=${model}
		# Copy generic files into directory and update
		outdir=$infDir/$pop/$model/$version/run_${SLURM_ARRAY_TASK_ID}
		mkdir -p $outdir

		cp $genericDir/$model.tpl $genericDir/$model.est $outdir # copy .est and .tpl files to outdir
		sed -i'' "s/SAMPLE_SIZE/$ss/g" $outdir/$model.tpl # sub in the sample size; note you need double quotes for variable to be expanded

		# Get sfs into inference directory and rename to match .est and .tpl files
		cp $sfsDir/${pop}_MAFpop0.obs $outdir/${header}_MAFpop0.obs # copy your sfs into the directory where you'll be doing the fsc inference
		cd $outdir
		fsc27093 -t ${header}.tpl -n 1000000 -m -e ${header}.est -M -L 60 -c${cores} -q

	done

done



