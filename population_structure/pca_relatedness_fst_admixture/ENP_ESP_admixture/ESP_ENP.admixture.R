# Title: Plot admixture output for all scenarios in all70 populations
# Author: Meixi Lin & Adjusted / Modified by Kaden Winspear @ CSUSM -> Eastern pacific Fin Whale project.
# May 2025 

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(ggplot2)
library(ggpubr)
library(dplyr)
# library(RColorBrewer)

# source functions
source('/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/plotting_config.R')

# def functions --------
read_admixture <- function(k, prefix, nrep, popmap, poporder) {
  # Q (the ancestry fractions) <-- What we need
  # P (the allele frequencies of the inferred ancestral populations).
  qfile = paste0(prefix, '.K', k,'.iter',1:nrep, '.Q')
  admaster = data.frame()
  for (ii in 1:length(qfile)) {
    if (!file.exists(qfile[ii])){
      print(qfile[ii])
      stop('File not exist')
    }
    ad = read.table(file = qfile[ii], stringsAsFactors = FALSE, header = FALSE)
    colnames(ad) = paste0("Cluster", 1:k)
    ad$run = ii
    ad = base::cbind(popmap, ad)
    admaster = base::rbind(admaster, ad)
  }
  admaster = admaster %>% dplyr::mutate(K = k)
  return(admaster)
}

format_admaster <- function(admaster, keeprun = c(2,4,6,8,10)) {
  K = unique(admaster$K)
  forplot = admaster %>%
    dplyr::filter(run %in% keeprun) %>%
    reshape2::melt(id.vars = c("SampleId", "PopId", "K", "run"), measure.vars = paste0('Cluster', 1:K)) %>%
    dplyr::mutate(K = paste('K =', K))
  return(forplot)
}

get_sampleorder <- function(popmap, poporder) {
  popmap$PopId = factor(popmap$PopId, levels = poporder)
  dt = popmap %>% dplyr::arrange(PopId)
  sampleorder = unique(dt$SampleId)
  return(sampleorder)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
mafcut = '05'
Klist = 1:6
poporder = c("ESP")

gdsfile = 'LDPruned_maf05_SNPs.ESPOnly_biallelic_all.gds'

workdir = paste('/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/ENP_ESP_admixture', sep = '/')
setwd(workdir)
indir = paste0('./Admixture_20250911/maf05/') # adjust path as needed
outdir = './derive_data/'
plotdir = './plots/'


dir.create(outdir, showWarnings = FALSE)
dir.create(plotdir, showWarnings = FALSE)

# load data --------

popmap = read.csv(file = "/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/ENP_ESP.popmap.csv", stringsAsFactors = F) # For all70 -> all70_fin_popmap.csv | LowGeno Data -> LowGenoRm_fin_popmap.csv | For DemoPass -> DemoPassInds_fin_popmap.csv


plinkname = read.table(file = '/data/shared/snigenda/finwhale_projects/fin_genomics/PopStructure/ENP_ESP_admixture/JointCalls_ESP_ENP_Only_biallelic_all_LDPruned_maf05_SA_mrf.nosex', header = FALSE, stringsAsFactors = FALSE)

if (!all(plinkname$V1 == popmap$SampleId)) {
  stop('Admixture file contains mismatching samples')
}

# load admixture dataframe ========
adprefix = paste0(indir, 'ESP_ENP_Only_biallelic_all_LDPruned')

dtlist = lapply(Klist, read_admixture, prefix = adprefix, nrep = 10, popmap = popmap, poporder = poporder)
forplotlist = lapply(dtlist, format_admaster, keeprun = 1:10)

# load CV data ========
adcv = read.csv(file = paste0(indir, 'Admixture_CV_ESP_ENP_Only_LLsummary_maf05.csv'), stringsAsFactors = FALSE)

# main --------
forplotQ = dplyr::bind_rows(forplotlist)
forplotQ$PopId = factor(forplotQ$PopId, levels = poporder)
forplotQ$SampleId = factor(forplotQ$SampleId, levels = get_sampleorder(popmap, poporder))

# full plot with all Ks and all 10 runs --------
pp <- ggplot(forplotQ, aes(x = SampleId, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  facet_grid(run ~ K) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(y = 'Ancestry Fraction') +
  theme_pubclean() +
  theme(
    axis.ticks.x = element_line(size = 0.3),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3), 
    axis.text.y = element_text(size = 7),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )


# Not combining plots since it was a pain ========
ggsave(filename = '10iterations_ENP_ESP_AdmixtureQ_all70_maf05.pdf', path = plotdir, plot = pp, height = 8, width = 12)


