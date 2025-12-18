# Title: Plot PCA by locations
# Author: Meixi Lin
# Date: Mon Jan 11 16:41:50 2021

# Modified / adjusted by Kaden Winspear @ CSUSM - > Eastern pacific fin whale project. 
# Date: April 2025


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)



library(dplyr)
library(ggplot2)
library(SNPRelate)
library(ggpubr)

# def functions --------
format_pca <- function(pca, popmap, poplevel) {
    pc.percent <- pca$varprop*100
    tab <- data.frame(sample.id = pca$sample.id,
                      pop = popmap[match(pca$sample.id, popmap$SampleId), poplevel],
                      PC1 = pca$eigenvect[,1], # the first eigenvector
                      PC2 = pca$eigenvect[,2], # the second eigenvector
                      PC3 = pca$eigenvect[,3], # the third eigenvector
                      PC4 = pca$eigenvect[,4], # the fourth eigenvector
                      stringsAsFactors = FALSE)
    lbls <- paste0("PC", 1:4, "(", format(pc.percent[1:4], digits=2), "%", ")")
    lims <- apply(pca$eigenvect[,1:4], 2, range)
    return(list(tab, lbls, lims))
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all70'
ref = 'Blue'
mafcut = '05'

workdir = paste('/Users/kadenwinspear/Documents/finwhale_proj/PopStructure', dataset, sep = '/')
outdir = './derive_data/'
plotdir = './plots/'

setwd(workdir)
sessionInfo()

dir.create(outdir)
dir.create(plotdir)
sessionInfo()
source('/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/plotting_config.R')

# load data --------
# final gdf file used
gdspath = "/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/Blue/data_files/JointCalls_all70_filterpass_biallelic_all_LDPruned_maf05.gds"
genofile <- snpgdsOpen(gdspath, readonly = TRUE)
popmap = read.csv(file = "/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/all70_fin_popmap.csv", stringsAsFactors = F)

sampleids = read.gdsn(index.gdsn(genofile, "sample.id"))
# check if the sampleids is the same order as popmap
if (!all(sampleids == popmap$SampleId)) {
    stop('GDS file contains mismatching samples')
}
# main --------
ENPids = sampleids[popmap$PopId %in% c('ENP', 'GOC', 'ESP')]

pcaout = snpgdsPCA(genofile, sample.id = ENPids, num.thread=8, autosome.only=FALSE, verbose = FALSE)
pcaformatted = format_pca(pcaout, popmap = popmap, poplevel = "PopId")
tab = pcaformatted[[1]];lbls = pcaformatted[[2]]; lims = pcaformatted[[3]]
tab$pop = factor(tab$pop, levels = poporder)

# plotting ========
pp1 <- ggplot(tab, aes(x=PC1, y =PC2, color= pop)) +
    geom_point(size = 2) +
    labs(x = lbls[1], y = lbls[2], color = 'Location') +
    coord_cartesian(xlim = lims[,1], ylim = lims[,2]) +
    scale_color_manual(values = mycolors) +
    theme_pubr() +
    theme(legend.position = "right")
legpp = as_ggplot(get_legend(pp1))
pp1 <- pp1 + theme(legend.position = "none")

pp2 <- ggplot(tab, aes(x=PC1, y =PC3, color= pop)) +
    geom_point(size = 2) +
    labs(x = lbls[1], y = lbls[3]) +
    coord_cartesian(xlim = lims[,1], ylim = lims[,3]) +
    scale_color_manual(values = mycolors) +
    theme_pubr() +
    theme(legend.position = "none")

pp3 <- ggplot(tab, aes(x=PC3, y =PC2, color= pop)) +
    geom_point(size = 2) +
    labs(x = lbls[3], y = lbls[2]) +
    coord_cartesian(xlim = lims[,3], ylim = lims[,2]) +
    scale_color_manual(values = mycolors) +
    theme_pubr() +
    theme(legend.position = "none")

pp <- ggarrange(pp1, pp3, pp2, legpp,
          labels = c("A", "C", "B"),
          ncol = 2, nrow = 2)

# output files --------
ggsave(filename = paste0("pca_maf", mafcut, "_all70_", today, ".pdf"), plot = pp, path = plotdir, height = 6, width = 6)

write.csv(tab, file = '/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/derive_data/PCA_results.csv')

# cleanup --------
# sink()
closefn.gds(genofile)
closeAllConnections()

