#!/bin/bash
#SBATCH --job-name=PopStructure      # Job name
#SBATCH --mem=20G                     # Total memory request (matches h_data=20G)
#SBATCH --time=23:00:00               # Wall time limit (HH:MM:SS)
#SBATCH --mail-type=ALL    # Notification events (abe = BEGIN,END,FAIL)
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --partition=CPU

# @usage        sbatch  step4_admixture_all50_20210316.sh
# Original Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin

# Adjusted by Kaden Winspear for pairwise comparisons between the populations. (ENP_GOC | ENP_ESP | ESP_GOC)

# Output:
# 1) 10 runs of admixture from K=1 to K=6
# 2) files with stats
# Notes:
# inputFile may be:
#      - a PLINK .bed file
#      - a PLINK "12" coded .ped file

############################################################
## import packages

module load bcftools/1.19
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate whalinGATK

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
DATASET='all70'
REF='Blue'
MAFCUT='05'
TODAY=$(date "+%Y%m%d")

HOMEDIR=/data/shared/snigenda/finwhale_projects/fin_genomics
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}
OUTDIR=${WORKDIR}/Admixture_All_Samples_Maf05
mkdir -p ${OUTDIR}

# Blue whale reference genome
REFERENCE=/data/shared/snigenda/finwhale_projects/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta

# sample to exclude
#EXCLUDE_SAMPLE=(GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")
#EXCLUDE_SAMPLE=(ESPCL01" "ESPCL02" "ESPCL03" "ESPCL04" "ESPCL05" "ESPCL06" "ESPCL07" "ESPCL08" "ESPCL09" "ESPCL10" "ESPCL11" "ESPCL12" "ESPCL13" "ESPCL14" "ESPCL15" "ESPCL16" "ESPCL17" "ESPCL18" "ESPCL19" "ESPCL20")
#EXCLUDE_SAMPLE=(ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15")
EXCLUDE_SAMPLE=()
# samples to include (not used)

#ENP_GOC_SAMPLES=ENPAK19,ENPAK20,ENPAK21,ENPAK22,ENPAK23,ENPAK24,ENPAK25,ENPAK26,ENPAK27,ENPAK28,ENPAK29,ENPAK30,ENPBC16,ENPBC17,ENPBC18,ENPCA01,ENPCA02,ENPCA03,ENPCA04,ENPCA05,ENPCA06,ENPCA07,ENPCA08,ENPCA09,ENPOR10,ENPOR11,ENPOR12,ENPOR13,ENPWA14,ENPWA15,GOC002,GOC006,GOC010,GOC025,GOC038,GOC050,GOC053,GOC063,GOC068,GOC071,GOC077,GOC080,GOC082,GOC086,GOC091,GOC100,GOC111,GOC112,GOC116,GOC125"
#ENP_ESP_SAMPLES=ENPAK19,ENPAK20,ENPAK21,ENPAK22,ENPAK23,ENPAK24,ENPAK25,ENPAK26,ENPAK27,ENPAK28,ENPAK29,ENPAK30,ENPBC16,ENPBC17,ENPBC18,ENPCA01,ENPCA02,ENPCA03,ENPCA04,ENPCA05,ENPCA06,ENPCA07,ENPCA08,ENPCA09,ENPOR10,ENPOR11,ENPOR12,ENPOR13,ENPWA14,ENPWA15,ESPCL01,ESPCL02,ESPCL03,ESPCL04,ESPCL05,ESPCL06,ESPCL07,ESPCL08,ESPCL09,ESPCL10,ESPCL11,ESPCL12,ESPCL13,ESPCL14,ESPCL15,ESPCL16,ESPCL17,ESPCL18,ESPCL19,ESPCL20"
ALL_SAMPLES="ENPAK19,ENPAK20,ENPAK21,ENPAK22,ENPAK23,ENPAK24,ENPAK25,ENPAK26,ENPAK27,ENPAK28,ENPAK29,ENPAK30,ENPBC16,ENPBC17,ENPBC18,ENPCA01,ENPCA02,ENPCA03,ENPCA04,ENPCA05,ENPCA06,ENPCA07,ENPCA08,ENPCA09,ENPOR10,ENPOR11,ENPOR12,ENPOR13,ENPWA14,ENPWA15,GOC002,GOC006,GOC010,GOC025,GOC038,GOC050,GOC053,GOC063,GOC068,GOC071,GOC077,GOC080,GOC082,GOC086,GOC091,GOC100,GOC111,GOC112,GOC116,GOC125,ESPCL01,ESPCL02,ESPCL03,ESPCL04,ESPCL05,ESPCL06,ESPCL07,ESPCL08,ESPCL09,ESPCL10,ESPCL11,ESPCL12,ESPCL13,ESPCL14,ESPCL15,ESPCL16,ESPCL17,ESPCL18,ESPCL19,ESPCL20"
# Original VCF file
INVCF="${WORKDIR}/JointCalls_all70_filterpass_bialleic_all_LDPruned_maf05_SA_mrF.vcf"
# Out prefix
OUTPREFIX="ALL_SAMPLES_LDPruned_maf05_SNPs"

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}; Maf cutoff = ${MAFCUT}; Input VCF = ${INVCF}"

############################################################
# start with step3 prune vcf to only include ENP samples
cd $OUTDIR

gatk --java-options "-Xmx4G" SelectVariants \
  -R $REFERENCE \
  $(for ss in "${EXCLUDE_SAMPLE[@]}"; do echo "--exclude-sample-name $ss"; done) \
  --exclude-filtered \
  --remove-unused-alternates \
  --exclude-non-variants \
  --restrict-alleles-to BIALLELIC \
  --select-type-to-include SNP \
  -V "${INVCF}" \
  -O "${OUTPREFIX}.vcf"


bcftools +counts ${OUTPREFIX}.vcf
conda deactivate

############################################################
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate fintools
############################################################
# convert to plink files
plink --vcf ${OUTPREFIX}.vcf --allow-extra-chr --recode 12 --out ${OUTPREFIX}

conda deactivate

############################################################
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate admixture
############################################################
# start to run admixture
# make a directory to store the P.Q files
mkdir -p Admixture_PQ

echo -e "K,iter,CVERROR,LL" > Admixture_CV_LLsummary_maf${MAFCUT}_${TODAY}.csv
for K in {1..6};do
    for i in {1..10};do
        # -s time the random seed to be generated from the current time
        # --cv In this default setting, the cross-validation procedure will perform 5-fold CV
        admixture --cv -s time -j8 ${OUTPREFIX}.ped ${K} | tee log_K${K}.iter${i}.out
        # need to move the files and keep all the runs
        mv ${OUTPREFIX}.${K}.Q Admixture_PQ/${OUTPREFIX}.K${K}.iter${i}.Q
        mv ${OUTPREFIX}.${K}.P Admixture_PQ/${OUTPREFIX}.K${K}.iter${i}.P
        # get the CV error and loglikelihood during each run
        CVERROR=$(awk '/^CV/ {print $4}' log_K${K}.iter${i}.out)
        LL=$(awk '/^Loglikelihood/ {print $2}' log_K${K}.iter${i}.out)
        echo -e "${K},${i},${CVERROR},${LL}" >> Admixture_CV_LLsummary_maf${MAFCUT}_${TODAY}.csv
    done
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"


conda deactivate
