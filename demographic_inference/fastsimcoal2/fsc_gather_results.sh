#!/bin/bash
#SBATCH --job-name=gatherResults
#SBATCH --output=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal/reports/fscresults.out.txt
#SBATCH --error=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal/reports/fscresults.err.txt
#SBATCH --time=00:09:00
#SBATCH --partition=HIMEM
#SBATCH --mem=8G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=winsp001@csusm.edu


# Gather up FSC results
# Author: Anabell Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Adjusted by Kaden Winspear @ CSUSM for slurm. Eastern pacific fin whale project.
# Usage: sbatch fsc_gather_results.sh


#Define directories ------------

wd=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal
infDir=$wd/3DModels


#Define Parameters -----------
models='AncestralSizeChange-Split-iso-assym'
pops="ENP-ESP"
version="Normal"


############################ MAIN ##############################

mkdir $infDir/resultsSummaries

for pop in $pops
do

	for model in $models
	do
		sumdir=$infDir/$pop/resultsSummaries/$model/$version
		mkdir -p $sumdir

		header=${model}

		outfile=$sumdir/${model}.results.Sumaries.csv

		header=`head -n1 $infDir/$pop/$model/$version/run_1/$model/*bestlhoods`
		echo -e "runNum\t$header" | tr "\\t" "," > $outfile

		for i in {1..100}
		do
			outdir=$infDir/$pop/$model/$version/run_${i}/$model 
			results=`grep -v [A-Z] $outdir/*bestlhoods`
			echo -e "$i\t$results"| tr "\\t" "," >> $outfile

		done
	done

done
