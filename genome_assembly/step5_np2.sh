#!/bin/bash
#SBATCH --job-name=np2_pipeline
#SBATCH --output=np2_pipeline_%j.out
#SBATCH --error=np2_pipeline_%j.err
#SBATCH --partition=HIMEM
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --mail-user=winsp001@csusm.com
#SBATCH --mail-type=ALL


#========WinBin===========
#        .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~
#=======Genomics===========

set -euo pipefail

module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate genome_ass

########## Directories ##########
ONT_READS="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/ONT_data/raw_data/finwhale_all.fastq.gz"
MEDAKA_ASM="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/medaka/consensus.fasta"
YAK_DB="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/yak/illumina.yak"

# Output dirs
MAPDIR="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/np2/map"
NP2DIR="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/np2"
QVDIR="/data/shared/snigenda/finwhale_projects/cetacean_genomes/gulf_genome_assembly/np2/qv"

mkdir -p "$MAPDIR" "$NP2DIR" "$QVDIR"

############################
# 5) Map ONT → Medaka asm #
############################
minimap2 -t 48 -x map-ont -a "$MEDAKA_ASM" "$ONT_READS" \
  | samtools sort -@ 16 -o "$MAPDIR/ont2medaka.bam"
samtools index "$MAPDIR/ont2medaka.bam"

#####################################
# 6) Illumina-guided polishing: NP2 #
#####################################
# Usage: nextPolish2 [OPTIONS] <HiFi.map.bam> <genome.fa[.gz]> <short.read.yak>...
nextPolish2 \
  -t 48 \
  -o "$NP2DIR/gulf.np2.fa" \
  "$MAPDIR/ont2medaka.bam" \
  "$MEDAKA_ASM" \
  "$YAK_DB"

##########################################
# 7) Re-measure QV on NP2-polished fasta #
##########################################
yak qv -t 48 -p -K3g -l100k "$YAK_DB" "$NP2DIR/gulf.np2.fa" > "$QVDIR/qv.np2.txt"

echo "Done."
echo "BAM:     $MAPDIR/ont2medaka.bam(.bai)"
echo "NP2 asm: $NP2DIR/gulf.np2.fa"
echo "QV:      $QVDIR/qv.np2.txt"
