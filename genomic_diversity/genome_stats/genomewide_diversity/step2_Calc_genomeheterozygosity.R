# Adjusted / modifed by Kaden Winspear @ CSUSM - > Eastern Pacific Fin Whale Project. 
# Date: April 4th 2025 

# Used in with wrapper. Counts heterozygous sites across genome.


rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(reshape2)

source("/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/plotting_config.R")

# def variables --------
dataset = "all70"
nsample = 70
ref = "Blue"
ncontig = 21 # 1 for Xchrom

today = format(Sys.Date(), "%Y%m%d")
datadate = "20250123" # the date the data generated

# go to that dataset
setwd('/data/shared/snigenda/finwhale_projects')
indir = './data/genome_stats/raw_data/'
outdir = '/data/shared/snigenda/finwhale_projects/fin_genomics/heterozygosity/genomewide_diversity/Xchromwide_estimates/'
dir.create(outdir)
sessionInfo()

# load data --------
#summaryfiles = paste0("/data/shared/snigenda/finwhale_projects/fin_genomics/heterozygosity/genomewide_diversity/count_sites_20250404/all70_Blue_", sprintf('%02d', seq(1,21)), "_sites_summary.txt") # For autosomes
summaryfiles = paste0("/data/shared/snigenda/finwhale_projects/fin_genomics/heterozygosity/genomewide_diversity/count_sites_20250404/all70_Blue_", sprintf('%02d', 22), "_sites_summary.txt") # For Xchrom

sumdt_list = lapply(summaryfiles, read.delim, stringsAsFactors = F)
sumdt = dplyr::bind_rows(sumdt_list)

blue_contiglist = read.csv(file = "/data/shared/snigenda/finwhale_projects/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/contig_list/blue_whale/blue_contig_summary.csv", row.names = 1, stringsAsFactors = F)

# main --------
sumdt = sumdt %>%
    dplyr::mutate(PopId = substr(SampleId, 1, 3),
                  SubPopId = ifelse(PopId == "ENP",
                                    substr(SampleId, 4, 5), PopId),
                  HomAltRealCount = HomAltCount - HomAltAllCount,
                  ContigHet = HetCount/CalledCount)

# Basic stats --------
PassCount = sum(unique(sumdt$PassCount))
# 890858824(all50_20201115)
TotalCount = sum(unique(sumdt$TotalCount))
# 2324429748(all50_20201115)
TotalRefLen = sum(blue_contiglist$LN)
# 2324429847
# Note the differences in TotalCount and TotalRefLen is mostly due to the errors during JointGenotyping such as too many alleles issues

# get the average heterozygosity
totaldt = sumdt %>%
    dplyr::group_by(SampleId, PopId) %>%
    dplyr::summarise(TotalHomRef = sum(HomRefCount),
                     TotalHomAlt = sum(HomAltCount),
                     TotalHomAltAll = sum(HomAltAllCount),
                     TotalHet = sum(HetCount),
                     TotalCalled = sum(CalledCount),
                     TotalMissing = sum(MissingCount),
                     TotalPass = sum(PassCount),
                     TotalgVCF = sum(as.numeric(TotalCount)),
                     .groups = 'drop') %>%
    dplyr::mutate(GenomeHet = TotalHet/TotalCalled) # genomewide heterozygosity

# output files --------
#write.csv(totaldt, file = paste0(outdir, "all70_genomewide_autosomal_heterozygosity_", today, ".csv"))
write.csv(totaldt, file = paste0(outdir, "all70_genomewide_XChrom_heterozygosity_", today, ".csv"))
#write.csv(sumdt, file = paste0(outdir, "all70_sites_summary_autosomal_combined_", today, ".csv"))
write.csv(sumdt, file = paste0(outdir, "all70_sites_summary_XChrom_combined_", today, ".csv"))

# cleanup --------


