#!/bin/bash
#SBATCH --job-name=subset_ROH
#SBATCH --time=16:00:00
#SBATCH --partition=HIMEM
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu

## This script is for subsetting bcftools roh output by individual, since the output files are very large
## Author: Sergio Nigenda & Kaden Winspear @ CSUSM -> Eastern pacific Fin whale project.
## The script also compresses the files, as they can be read into R when gzipped


cd /data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/ESP

file=ESP_concat_fwhale_roh_bcftools_G30_ACANGT.out

tail -n +6 ${file} > tmp.txt


awk '$2 == "ESPCL01" { print }' tmp.txt > ESPCL01_${file}
awk '$2 == "ESPCL02" { print }' tmp.txt > ESPCL02_${file}
awk '$2 == "ESPCL03" { print }' tmp.txt > ESPCL03_${file}
awk '$2 == "ESPCL04" { print }' tmp.txt > ESPCL04_${file}
awk '$2 == "ESPCL05" { print }' tmp.txt > ESPCL05_${file}
awk '$2 == "ESPCL06" { print }' tmp.txt > ESPCL06_${file}
awk '$2 == "ESPCL07" { print }' tmp.txt > ESPCL07_${file}
awk '$2 == "ESPCL08" { print }' tmp.txt > ESPCL08_${file}
awk '$2 == "ESPCL09" { print }' tmp.txt > ESPCL09_${file}
awk '$2 == "ESPCL10" { print }' tmp.txt > ESPCL10_${file}
awk '$2 == "ESPCL11" { print }' tmp.txt > ESPCL11_${file}
awk '$2 == "ESPCL12" { print }' tmp.txt > ESPCL12_${file}
awk '$2 == "ESPCL13" { print }' tmp.txt > ESPCL13_${file}
awk '$2 == "ESPCL14" { print }' tmp.txt > ESPCL14_${file}
awk '$2 == "ESPCL15" { print }' tmp.txt > ESPCL15_${file}
awk '$2 == "ESPCL16" { print }' tmp.txt > ESPCL16_${file}
awk '$2 == "ESPCL17" { print }' tmp.txt > ESPCL17_${file}
awk '$2 == "ESPCL18" { print }' tmp.txt > ESPCL18_${file}
awk '$2 == "ESPCL19" { print }' tmp.txt > ESPCL19_${file}
awk '$2 == "ESPCL20" { print }' tmp.txt > ESPCL20_${file}

wait

gzip ESP*
# gzip MN*
gzip ${file}
rm tmp.txt
