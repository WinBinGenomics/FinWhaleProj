#! /bin/bash
#SBATCH --job-name=exon_dist
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --partition=CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=winsp001@csusm.edu

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools

HOMEDIR=/data/shared/snigenda/finwhale_projects
WORKDIR=${HOMEDIR}/fin_genomics/demographic_inference/Neutral_regions

REFDIR=${HOMEDIR}/cetacean_genomes/blue_whale_genome/
GFF=${REFDIR}/annotation/GCF_009873245.2_mBalMus1.pri.v3_genomic.gff
CDS_REGIONS=${REFDIR}/annotation/GCF_009873245.2_mBalMus1.pri.v3_genomic_cds.bed
EXONS=${REFDIR}/annotation/GCF_009873245.2_mBalMus1.pri.v3_genomic_exons.bed
CDRS=${REFDIR}/annotation/GCF_009873245.2_mBalMus1.pri.v3_allCodingRegions_merged.bed
hqSites=${WORKDIR}/HiQualCoords/all_HQCoords_sorted_merged.bed

mkdir -p ${WORKDIR}/DistanceFromExons

############################

### Get exonic regions or coding regions (do only once) and make sure is sorted
#awk '$3=="exon"' GCF_009873245.2_mBalMus1.pri.v3_genomic.gff | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > GCF_009873245.2_mBalMus1.pri.v3_genomic_exons.bed
#awk '$3=="CDS"' GCF_009873245.2_mBalMus1.pri.v3_genomic.gff | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > GCF_009873245.2_mBalMus1.pri.v3_genomic_cds.bed

### Then concatenate the .bed files & sort so all exons and cds regions are together.
#cat GCF_009873245.2_mBalMus1.pri.v3_genomic_cds.bed GCF_009873245.2_mBalMus1.pri.v3_genomic_exons.bed  | sort -k1,1 -k2,2n   >  GCF_009873245.2_mBalMus1.pri.v3_allCodingRegions.bed

### finally use bedtools merge to be left with a completely merged file consisting of non overlapping exon and CDS regions.
#bedtools merge -i GCF_009873245.2_mBalMus1.pri.v3_allCodingRegions.bed > GCF_009873245.2_mBalMus1.pri.v3_allCodingRegions_merged.bed

#############################

## There are discrepancies between the two methods, Annabel's approach results in 499895 lines, while mine has 499380 lines.
## I ended up using mine because it seems more specific to the column were the features are defined in gff files.

bedtools closest -d -a ${hqSites} -b ${EXONS} > ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt
bedtools closest -d -a ${hqSites} -b ${CDRS} > ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt


# Exploring and choosing different distances

# The last column (8) is the distance; I want it to be at least 10,000, and want to keep track of the distance. Collect all that are >10,000 away.
# Pick the ones with specific distance (awk) from exons, we tried 3 different distances: 10Kb, 20Kb, 50Kb.
awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromExons.0based.bed
awk -F'\t' '{OFS="\t";if($8>20000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromExons.0based.bed
awk -F'\t' '{OFS="\t";if($8>50000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min50kb_DistFromExons.0based.bed
# Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)

awk -F'\t' '{OFS="\t";if($7>10000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromCDR.0based.bed
awk -F'\t' '{OFS="\t";if($7>20000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromCDR.0based.bed
awk -F'\t' '{OFS="\t";if($7>50000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min50kb_DistFromCDR.0based.bed

> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
for i in `ls ${WORKDIR}/DistanceFromExons/*DistFromExons.0based.bed`
do
echo $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
done

> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
for i in `ls ${WORKDIR}/DistanceFromExons/*DistFromCDR.0based.bed`
do
echo $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
done

conda deactivate

