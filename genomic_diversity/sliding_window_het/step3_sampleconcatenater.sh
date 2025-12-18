#!/bin/bash

# author -> Kaden winspear @ CSUSM: Eastern Pacific Fin whale project
# Date: April 2025

# This script will be used to cocatenate the create sample specific files. Each file should have 21 (or the amount of VCfs) files. This script ensures that each sample is combined so
# we are left with 70 (or the amount of samples) files containing window het data.

# Now that we have a file for each chrom of each sample we need to concatenate them for plotting / calculataion purposes. 

# Output: One comnined file for each sample. So 70 total, all containing the  21 autosomes for that respective individual.

cd /data/shared/snigenda/finwhale_projects/fin_genomics/diversity/sliding_window/1mbwins/

# Sample Ids

samples=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "ESPCL01" "ESPCL02" "ESPCL03" "ESPCL04" "ESPCL05" "ESPCL06" "ESPCL07" "ESPCL08" "ESPCL09" "ESPCL10" "ESPCL11" "ESPCL12" "ESPCL13" "ESPCL14" "ESPCL15" "ESPCL16" "ESPCL17" "ESPCL18" "ESPCL19" "ESPCL20" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")


# Loop over each sample
for sample in "${samples[@]}"; do
    # Initialize an empty output file for each sample
    output_file="combined_${sample}.txt"

    # Loop through the files for each sample and concatenate them
    first=true
    for file in ${sample}*.txt; do
        if [ "$first" = true ]; then
            # If it's the first file, include the header and concatenate
            head -n 1 "$file" > "$output_file"
            tail -n +2 "$file" >> "$output_file"
            first=false
        else
            # For subsequent files, skip the header and append the data
            tail -n +2 "$file" >> "$output_file"
        fi
    done

    echo "Finished combining files for $sample. Output saved in $output_file"

done

# clean up the dir

rm *temp.txt*
