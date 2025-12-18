# Title: Plot the ape neighbor joining tree
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin
# Date: Mon Mar  1 02:15:28 2021
# Example: Rscript --vanilla <homedir>/fin_whale/scripts_analyses/PopStructure/f50b4/step3.2_ApePhylogeny_f50b4_plot_20210301.R '<homedir>/finwhale/analyses/PopStructure/f50b4/Minke' 'f50b4_pass_bialleic_all_LDPruned_maf05'

# IMPORTANT NOTE: The input rows (sample names) for ape::nj need to be named with values other than 1,2,3... The dist.gene function will confuse the dimmension of genotype
# IMPORTANT NOTE: Here the input matrix should be on the dosage of alternative alleles (when you use the dosage of reference allele, tree looks differently.)
# IMPORTANT NOTE: The ggtree is pretty outdated in bioconda. https://github.com/YuLab-SMU/ggtree/issues/91

# Adjusted by Kaden Winspear @ CSUSM -> Eastern Fin Whale Pacific project. 
# No outgroup in this version. 

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)

source('/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/plotting_config.R')


# def functions --------
# get maf
get_maf <- function(outprefix) {
  maf = stringr::str_split(outprefix, pattern = '_')[[1]]
  maf = maf[length(maf)]
  return(maf)
}

# empty legend
theme_blank_legend <- theme(
  legend.position = 'bottom',
  panel.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent", colour = NA), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent", colour = NA), # get rid of legend panel bg
  legend.key = element_rect(fill = "transparent", colour = NA) # get rid of key legend fill, and of the surrounding
)

plot_trees <- function(tree_nj, popmap, main, nodecutoff = 750) {
  metadt0 <- dplyr::left_join(data.frame('label' =  tree_nj$tip.label), popmap,
                              by = c('label' = 'SampleId'))
  p1_rt <- ggtree(tree_nj)
  p1dt <- p1_rt$data
  metadt <- dplyr::full_join(p1dt, metadt0, by = 'label') %>%
    dplyr::mutate(boot_strap_support = ifelse(isTip == FALSE & label != 'Root', label, NA)) %>%
    dplyr::mutate(nodebt = ifelse(as.integer(boot_strap_support) >= nodecutoff, boot_strap_support, NA))
  p1_rt$data <- metadt
  outpp <- p1_rt +
    geom_tiplab(aes(color = PopId), show.legend = FALSE) +
    geom_point(aes(size = as.integer(nodebt), fill = as.integer(nodebt)), shape = 21, alpha = 0.8) +
    scale_size_binned(range = c(1, 8)) +
    scale_fill_gradient(low = 'white', high = 'red') +
    scale_color_manual(values = c(speccolors[-5], loccolors)) +
    labs(title = main, size = 'Bootstrap Support', fill = '') +
    coord_cartesian(xlim = c(min(metadt$x), 1.1 * max(metadt$x))) +
    theme_blank_legend
  outlist = list(tree_nj, metadt, outpp)
  return(outlist)
}

# plot the density trees
plot_densitrees <- function(btrees, alpha = 0.1, main) {
  pp <- ggdensitree(btrees, alpha = alpha, colour = 'steelblue') +
    geom_tiplab(size = 3) +
    labs(title = main)
  return(pp)
}

# def variables --------

workdir <- '/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/Blue'
outprefix <- 'lowGenoRemove_phylo'

today = format(Sys.Date(), "%Y%m%d")
datadate = '20250411'
nbt = 1000
maf = get_maf(outprefix)
outdir = './derive_phylogeny/'
plotdir = './plots/'

setwd(workdir)
sessionInfo()
# dir.create(outdir)
dir.create(plotdir)

# load data --------
treespair <- readRDS(file = paste0(outdir, outprefix, '_treespair_20250411.rds'))
treesperc <- readRDS(file = paste0(outdir, outprefix, '_treesperc_20250411.rds'))

# get the population names
popmap = read.csv(file = '/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/all70_fin_popmap.csv', stringsAsFactor = FALSE)

# main --------

# plot the ape trees ========
pp_pair <- plot_trees(tree_nj = treespair[[1]], popmap, main = paste0('NJ pairwise ', maf))
pp_perc <- plot_trees(tree_nj = treesperc[[1]], popmap, main = paste0('NJ percentage ', maf))

ggsave(filename = paste0(outprefix, '_UttreePair_LowGenoRemove', today, '.pdf'), path = plotdir, plot = pp_pair[[3]], width = 8, height = 8)
ggsave(filename = paste0(outprefix, '_UttreePerc_LowGenoRemove', today, '.pdf'), path = plotdir, plot = pp_perc[[3]], width = 8, height = 8)

# Modification: Add source data ========
# Date: Mon Jan 16 12:02:26 2023
pp_pairlist <- plot_trees(tree_nj = treespair[[1]], popmap, main = paste0('NJ pairwise ', maf))
mytree <- pp_pairlist[[1]]
outdt <- pp_pairlist[[2]]

ape::write.tree(mytree, file = '/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/derive_data/FigS5_lowgenoremove-tree.newick.txt')
write.csv(outdt, file = '/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/derive_data/tree_lowgenoremove.newick.csv')

# cleanup --------
date()
closeAllConnections()
