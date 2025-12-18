#!/bin/bash
#SBATCH --job-name=pi
#SBATCH --output=pi.out
#SBATCH --time=00:05:00
#SBATCH --mem=1G

#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Snigenda workshop
# Author:Kaden Winspear - Dec. 2025
# Usage:sbatch step2_average_pi_calc.sh
# Purpose: This script takes in Pixy output & calculates average pi per pop.

workdir="/data/shared/snigenda/finwhale_projects/filteredvcf/workshop/Diversity_analysis/pixy"
input="${workdir}/Chrom19thru21_pi.txt"

cd ${workdir}

# extract header and add to separated txt files
head -n 1 ${input} > ENP_pi.txt
head -n 1 ${input} > GOC_pi.txt

# split by population (column 1)
awk 'NR>1 && $1=="ENP"' ${input} >> ENP_pi.txt
awk 'NR>1 && $1=="GOC"' ${input} >> GOC_pi.txt

# compute means of avg_pi (column 5)
awk 'NR>1 {sum+=$5;n++} END {print "ENP mean pi:",sum/n}' ENP_pi.txt > mean_pi.txt
awk 'NR>1 {sum+=$5;n++} END {print "GOC mean pi:",sum/n}' GOC_pi.txt >> mean_pi.txt
