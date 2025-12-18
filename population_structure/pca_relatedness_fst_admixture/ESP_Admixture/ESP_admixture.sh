#!/bin/bash
#SBATCH --mem=20G
#SBATCH --time=23:00:00
#SBATCH --partition=HIMEM
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=winsp001@csusm.edu
#SBATCH --mail-type=ALL

# @version      v0
# @usage        qsub step4_admixture_all50_20210316.sh
# @description  Performs admixture analyses on the LDpruned sites (MAF cutoff = NA/0.05/0.10)
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin
# Date: Thu Mar 18 15:22:26 2021
# Output:
# 1) 10 runs of admixture from K=2 to K=6
# 2) files with stats
# Notes:
# inputFile may be:
#      - a PLINK .bed file
#      - a PLINK 12 coded .ped file

############################################################
## import packages
set -eo pipefail

############################################################
## def environments

module load R/4.3.1+Bioconductor
module load anaconda3/2023.09
source /cm/shared/apps/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate admixture

############################################################
## def variables
DATASET='ESPOnly'
MAFLISTS=('05')
TODAY=$(date "+%Y%m%d")

HOMEDIR="/data/shared/snigenda/finwhale_projects"
WORKDIR="${HOMEDIR}/fin_genomics/PopStructure/ESP_admixture"

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID}."
for MAFCUT in ${MAFLISTS[@]}; do
    echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Start."
    OUTDIR=${WORKDIR}/Admixture_${TODAY}/maf${MAFCUT}
    mkdir -p ${OUTDIR}
    cd ${OUTDIR}
    # filenames
    OUT=JointCalls_ESPOnly_biallelic_all_LDPruned_maf05_SA_mrf
    SHORTOUT=ESPOnly_biallelic_all_LDPruned
    PED=${WORKDIR}/${OUT}.ped
    echo -e "K,iter,CVERROR,LL" > Admixture_CV_DemoPassInds_LLsummary_maf05.csv
    for K in {1..6};do
        for i in {1..10};do
            # -s time the random seed to be generated from the current time
            # --cv In this default setting, the cross-validation procedure will perform 5-fold CV
            admixture --cv -s time -j8 ${PED} ${K} | tee log_K${K}.iter${i}.out
            mv ${OUT}.${K}.Q ${SHORTOUT}.K${K}.iter${i}.Q
            mv ${OUT}.${K}.P ${SHORTOUT}.K${K}.iter${i}.P
            # get the CV error and loglikelihood during each run
            CVERROR=$(awk '/^CV/ {print $4}' log_K${K}.iter${i}.out)
            LL=$(awk '/^Loglikelihood/ {print $2}' log_K${K}.iter${i}.out)
            echo -e "${K},${i},${CVERROR},${LL}" >> Admixture_CV_ESPOnly_LLsummary_maf${MAFCUT}.csv
        done
    done
    echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = 05; Done."
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOB_ID} Done"

