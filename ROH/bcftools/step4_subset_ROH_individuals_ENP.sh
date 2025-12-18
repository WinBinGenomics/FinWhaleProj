#!/bin/bash
#SBATCH --job-name=subset_ROH
#SBATCH --time=16:00:00
#SBATCH --partition=HIMEM
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=winsp001@csusm.edu


## This script is for subsetting bcftools roh output by individual, since the output files are very large.
## Author: Sergio Nigenda & Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale project.
## The script also compresses the files, as they can be read into R when gzipped

cd /data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/ENP


file=ENP_concat_fwhale_roh_bcftools_G30_ACANGT.out

tail -n +6 ${file} > tmp.txt


awk '$2 == "ENPAK19" { print }' tmp.txt > ENPAK19_${file}
awk '$2 == "ENPAK20" { print }' tmp.txt > ENPAK20_${file}
awk '$2 == "ENPAK21" { print }' tmp.txt > ENPAK21_${file}
awk '$2 == "ENPAK22" { print }' tmp.txt > ENPAK22_${file}
awk '$2 == "ENPAK23" { print }' tmp.txt > ENPAK23_${file}
awk '$2 == "ENPAK24" { print }' tmp.txt > ENPAK24_${file}
awk '$2 == "ENPAK25" { print }' tmp.txt > ENPAK25_${file}
awk '$2 == "ENPAK26" { print }' tmp.txt > ENPAK26_${file}
awk '$2 == "ENPAK27" { print }' tmp.txt > ENPAK27_${file}
awk '$2 == "ENPAK28" { print }' tmp.txt > ENPAK28_${file}
awk '$2 == "ENPAK29" { print }' tmp.txt > ENPAK29_${file}
awk '$2 == "ENPAK30" { print }' tmp.txt > ENPAK30_${file}
awk '$2 == "ENPBC16" { print }' tmp.txt > ENPBC16_${file}
awk '$2 == "ENPBC17" { print }' tmp.txt > ENPBC17_${file}
awk '$2 == "ENPBC18" { print }' tmp.txt > ENPBC18_${file}
awk '$2 == "ENPCA02" { print }' tmp.txt > ENPCA02_${file}
awk '$2 == "ENPCA03" { print }' tmp.txt > ENPCA03_${file}
awk '$2 == "ENPCA04" { print }' tmp.txt > ENPCA04_${file}
awk '$2 == "ENPCA05" { print }' tmp.txt > ENPCA05_${file}
awk '$2 == "ENPCA06" { print }' tmp.txt > ENPCA06_${file}
awk '$2 == "ENPCA07" { print }' tmp.txt > ENPCA07_${file}
awk '$2 == "ENPCA08" { print }' tmp.txt > ENPCA08_${file}
awk '$2 == "ENPOR10" { print }' tmp.txt > ENPOR10_${file}
awk '$2 == "ENPOR11" { print }' tmp.txt > ENPOR11_${file}
awk '$2 == "ENPOR13" { print }' tmp.txt > ENPOR13_${file}
awk '$2 == "ENPWA14" { print }' tmp.txt > ENPWA14_${file}
awk '$2 == "ENPWA15" { print }' tmp.txt > ENPWA15_${file}


wait

gzip ENP*
# gzip MN*
gzip ${file}
rm tmp.txt
