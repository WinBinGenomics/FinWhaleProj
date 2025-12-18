# Author: Meixi Lin
# Date: Sat Jan  1 16:41:41 2022

#Adjusted / Modified by Kaden Winspear @ CSUSM - > Eastern Pacific Fin Whale project. 
# Modification: Adjusted for ENP/GOC/ESP only. No subpopulations. added ESP analysis. 
# Date: April 2025.

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(dunn.test)

source('/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/plotting_config.R')

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

setwd('/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70')
outdir = './derive_data/revisions_Kinship/'
plotdir = './plots/revisions_Kinship/'

dir.create(outdir, showWarnings = FALSE)
dir.create(plotdir, showWarnings = FALSE)

dataset = 'all70'
ref = 'Blue'
mafcut = '05' # mafcutoff across the population used in the input gdsfile (NOT the maf setting used in the other samples)
missing = 0.2
gdsfile = './data_files/JointCalls_all70_filterpass_biallelic_all_LDPruned_maf05.gds'

sessionInfo()
getwd()

# color palettes ========
pop_colors = loccolors[c('ENP','GOC','ESP')] # Original color selection

# def function --------
estimate_kinship <- function(genofile, targetpop) {
  sampleids = popmap[popmap[,"PopId"] == targetpop,'SampleId']
  # Estimating IBD Using PLINK method of moments (MoM)
  # Missing rate and MAF set as the same in the tutorial (http://corearray.sourceforge.net/tutorials/SNPRelate/#f_st-estimation) and same as Paulina's previous settings
  print(paste0('Target population = ', targetpop))
  ibddt <- snpgdsIBDMoM(genofile, sample.id = sampleids, verbose = TRUE,
                        maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE)
  # k0        the probability of sharing ZERO alleles
  # k1        the probability of sharing ONE alleles
  # kinship   kinship coefficient
  ibddt.coeff <- snpgdsIBDSelection(ibddt)
  ibddt.coeff = dplyr::left_join(ibddt.coeff, popmap, by = c('ID1' = 'SampleId'))
  return(ibddt.coeff)
}

plot_kinship <- function(ibd.coeff) {
  # prepare data
  ibd.coeff$ID1 = factor(ibd.coeff$ID1, levels = popmap$SampleId)
  ibd.coeff$ID2 = factor(ibd.coeff$ID2, levels = popmap$SampleId)
  ibd.coeff$PopId = factor(ibd.coeff$PopId, levels = c('ENP','GOC','ESP'))
  pp1 <- ggplot(ibd.coeff, aes(x=ID1, y=ID2, fill=kinship)) +
    geom_tile() + 
    theme_light() +
    scale_fill_gradient(low = '#1F618D', high = '#2ECC71') +
    theme(axis.text.x=element_text(angle=90),
          aspect.ratio = 1)
  
  # get a boxplot
  pp2 <- ggplot(ibd.coeff, aes(x = PopId, y = kinship, color = PopId)) +
    geom_boxplot() +
    scale_color_manual(values = pop_colors) +
    theme_light() +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          aspect.ratio = 1)
  
  # set for inset (original positioning kept)
  pp <- pp1 +
    annotation_custom(
      ggplotGrob(pp2),
      xmin = 24, xmax = 44, ymin = 1, ymax = 21
    )
  return(pp)
}

kinship_stats <- function(ibd) {
  print('#######################################################################')
  # excluding zeros
  ibdsum = ibd %>%
    dplyr::group_by(PopId) %>%
    dplyr::summarise(mean_kinship = mean(kinship),
                     mean_kinship_plus0 = sum(kinship)/sum(kinship>0),
                     .groups = 'drop')
  print(ibdsum)
  print(ibd[,c('PopId','kinship')] %>% split(., ibd[,'PopId']) %>% purrr::map(summary))
  
  # run dunn tests and kruskal tests
  ibd$PopId = factor(ibd$PopId, levels = c('ENP','GOC','ESP'))
  print(kruskal.test(ibd[,'kinship'] ~ ibd[,'PopId']))
  if (nlevels(ibd$PopId) > 2) {
    dunn.test::dunn.test(x = ibd[,'kinship'], g = ibd[,'PopId'], method = 'none')
    dunn.test::dunn.test(x = ibd[,'kinship'], g = ibd[,'PopId'], method = 'bonferroni')
    dunn.test::dunn.test(x = ibd[,'kinship'], g = ibd[,'PopId'], method = 'bh')
  }
  return(ibdsum)
}

# load data --------
genofile = SNPRelate::snpgdsOpen(filename = paste0(ref, '/', gdsfile))

#List of sample ids
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# sample designations (modified for target populations)
popmap = read.csv(file = '/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/all70_fin_popmap.csv', stringsAsFactors = FALSE) %>%
  dplyr::filter(PopId %in% c('ENP','GOC','ESP')) %>% # Filter to target populations
  dplyr::arrange(PopId) # Original ordering preserved

# main --------
# estimate IBD within populations ========
pop_ibdlist <- lapply(c('ENP','GOC','ESP'), function(xx) {
  ibd <- estimate_kinship(genofile, targetpop = xx)
})

pop_ibd <- dplyr::bind_rows(pop_ibdlist)

pp <- plot_kinship(pop_ibd)
ggsave(filename = paste0(plotdir, 'kinship_populations_', today, '.pdf'), plot = pp, height = 8, width = 8)

# get summary ########
pop_mean <- kinship_stats(ibd = pop_ibd)

# cleanup --------
closefn.gds(genofile)
date()

# save image --------
save.image(file = paste0(outdir, 'revisions_Relatedness_all70_bygroup_', today, '.RData'))
closeAllConnections()

# load the image to output the data --------
load('./derive_data/revisions_Kinship/revisions_Relatedness_all70_bygroup_20250410.RData')

# Final output handling (original structure preserved)
outplotdata = pop_ibd %>%
  dplyr::select(ID1, ID2, kinship, PopId)

write.csv(outplotdata,file = '/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/derive_data/revisions_Kinship/Kinship_data.csv')
