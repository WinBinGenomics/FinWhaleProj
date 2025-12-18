#!/bin/bash

# author -> Kaden Winspear @ CSUSM Eastern Pacific Fin Whale project.
# Date April 2025.

# purpose: Following the python script you are left with 21 files consisiting of the window data for all 70 samples. For plotting purposes we must parse these 21 files, extracting
# individuals and making new files that are sample specific. Following this script you were hav 21 sample specific files (one for each chromosome) which will then be concatenated using
# the next script. Sample_concatenater.

scriptdir='/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/genomic_diversity/sliding_window_het/'

cd  /data/shared/snigenda/finwhale_projects/fin_genomics/diversity/sliding_window/1mbwins/

# loops through each file matching the file name.
for file in JointCalls_08_b_VariantFiltration_*_passing.vcf.gz_het_1000000win_1000000step.txt; do
    echo "Processing $file"

    file_number=$(echo "$file" | sed -E 's/.*VariantFiltration_([0-9]+)_passing.*/\1/') #extracts chromosome number from the filename w/ sed. sed uses the search / replace to identify chrom ID.

    # runs parser script on files.
    "${scriptdir}/step2_parser.awk" "$file"

    echo "Finished processing $file"
done
