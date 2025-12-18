#!/bin/bash

# Author: Kaden Winspear @ CSUSM - > Eastern pacific fin whale project.
# Usage: sbatch pixy_calc.sh -> used to calculate / average out genomewide pi.

workdir='/data/shared/snigenda/finwhale_projects/fin_genomics/heterozygosity/genomewide_diversity/pixy_estimates/genomewide/'

cd ${workdir}/pi_estimates


# Step 1: Collect data (skip headers, keep population and pi columns)
for file in Whole_chrom_Analysis_*_pi.txt; do
  tail -n +2 "$file" | cut -f1,5 >> temp_data.txt
done

# Step 2: Get unique population names (e.g., ENP, ESP, GOC)
pops=$(cut -f1 temp_data.txt | sort -u)

# Step 3: Calculate and print averages
echo "pop    avg_pi"  # Header
for pop in $pops; do
  # Extract pi values for this population
  pi_values=$(grep "^$pop" temp_data.txt | cut -f2)

  # Sum the values
  sum=$(echo "$pi_values" | paste -sd+ | bc)

  # Count entries
  count=$(echo "$pi_values" | wc -l)

  # Calculate average (with 9 decimal places)
  average=$(echo "scale=9; $sum / $count" | bc)

  # Print result
  echo "$pop    $average" >> wholegenome_pi.txt
done

# Cleanup
rm temp_data.txt
