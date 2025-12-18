#!/bin/bash
#SBATCH --job-name=gatherResults
#SBATCH --output=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal/reports/fscresults.out.txt
#SBATCH --error=/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/models/fastsimcoal/reports/fscresults.err.txt
#SBATCH --time=00:09:00
#SBATCH --partition=CPU
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

# Define models
models='Split-symm-m01-asym-m12'

############################ MAIN ##############################

mkdir -p "$infDir/resultsSummaries"

for model in $models
do
    # Summary directory
    sumdir="$infDir/resultsSummaries/$model"
    mkdir -p "$sumdir"

    outfile="$sumdir/${model}.results.Summaries.csv"

    # Grab header from first run
    header_line=$(head -n1 "$infDir/$model/run_1/$model"/*bestlhoods)

    # Write CSV header
    echo -e "runNum\t$header_line" | tr "\\t" "," > "$outfile"

    # Loop through runs
    for i in {1..100}
    do
        outdir="$infDir/$model/run_${i}/$model"

        results=$(grep -v '[A-Z]' "$outdir"/*bestlhoods)

        echo -e "$i\t$results" | tr "\\t" "," >> "$outfile"
    done
done
