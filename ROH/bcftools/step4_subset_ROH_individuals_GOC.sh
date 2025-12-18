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


cd /data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/GOC

file=GOC_concat_fwhale_roh_bcftools_G30_ACANGT.out

tail -n +6 ${file} > tmp.txt


awk '$2 == "GOC002" { print }' tmp.txt > GOC002_${file}
awk '$2 == "GOC006" { print }' tmp.txt > GOC006_${file}
awk '$2 == "GOC025" { print }' tmp.txt > GOC025_${file}
awk '$2 == "GOC038" { print }' tmp.txt > GOC038_${file}
awk '$2 == "GOC050" { print }' tmp.txt > GOC050_${file}
awk '$2 == "GOC053" { print }' tmp.txt > GOC053_${file}
awk '$2 == "GOC063" { print }' tmp.txt > GOC063_${file}
awk '$2 == "GOC068" { print }' tmp.txt > GOC068_${file}
awk '$2 == "GOC071" { print }' tmp.txt > GOC071_${file}
awk '$2 == "GOC080" { print }' tmp.txt > GOC080_${file}
awk '$2 == "GOC082" { print }' tmp.txt > GOC082_${file}
awk '$2 == "GOC086" { print }' tmp.txt > GOC086_${file}
awk '$2 == "GOC091" { print }' tmp.txt > GOC091_${file}
awk '$2 == "GOC100" { print }' tmp.txt > GOC100_${file}
awk '$2 == "GOC111" { print }' tmp.txt > GOC111_${file}
awk '$2 == "GOC112" { print }' tmp.txt > GOC112_${file}
awk '$2 == "GOC116" { print }' tmp.txt > GOC116_${file}
awk '$2 == "GOC125" { print }' tmp.txt > GOC125_${file}


wait

gzip GOC*
# gzip MN*
gzip ${file}
rm tmp.txt
