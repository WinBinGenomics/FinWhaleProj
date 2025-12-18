# Title: Plot ZooRoH output along chromosomes
# Author: Meixi Lin (modified by Kaden)
# Date: Mon May 28 2025

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH")

library(dplyr)
library(ggplot2)
library(ggpubr)

# mapping scaffold IDs to chromosome numbers --------
chromosome_map <- c(
  "NC_045785.1" = "1", "NC_045786.1" = "2", "NC_045787.1" = "3",
  "NC_045788.1" = "4", "NC_045789.1" = "5", "NC_045790.1" = "6",
  "NC_045791.1" = "7", "NC_045792.1" = "8", "NC_045793.1" = "9",
  "NC_045794.1" = "10", "NC_045795.1" = "11", "NC_045796.1" = "12",
  "NC_045797.1" = "13", "NC_045798.1" = "14", "NC_045799.1" = "15",
  "NC_045800.1" = "16", "NC_045801.1" = "17", "NC_045802.1" = "18",
  "NC_045803.1" = "19", "NC_045804.1" = "20", "NC_045805.1" = "21"
)

# def functions --------
load_contiglist <- function(contigfile) {
  contiglist = read.csv(file = contigfile, row.names = 1, stringsAsFactors = FALSE)
  contiglist = contiglist %>%
    dplyr::mutate(gnPOS = genomewide_coord - 1,
                  chrom = chromosome_map[as.character(SN)]) %>%
    dplyr::select(chrom, gnPOS)
  colnames(contiglist) = c('chrom', 'gnPOS')
  return(contiglist)
}

# def variables --------
indir = "./derived_data/"
outdir = "./derived_data/"
plotdir = './plots/'
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

today = format(Sys.Date(), "%Y%m%d")
# total autosomal length
genomelen = 2239549461

# load data --------
zooroh = readRDS(file = './derived_data/ROH_zooroh_summary_final_20250602.rds')
contiglist = load_contiglist(
  contigfile = "/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/blue_contig_summary_with_coord.csv"
)

# prepare axis breaks and labels --------
chrom_breaks = contiglist$gnPOS
chrom_labels = contiglist$chrom

# ROH categories --------
rohlens = c(0.1, 1, 5) * 1e+6
rohcats = c('0.1_1', '1_5', '5_Inf')
lenslab = c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb')

sessionInfo()

# main plotting --------
# set colors by category
rohcatbrewer = RColorBrewer::brewer.pal(name = "Reds", n = 9)[c(5,7,9)]
names(rohcatbrewer) = unique(zooroh$rohcat)

# panel 1: all ROH segments
pp2 <- ggplot(zooroh) +
  geom_segment(aes(x = gnstart, xend = gnend,
                   y = sample, yend = sample,
                   color = rohcat), size = 3) +
  scale_color_manual(values = rohcatbrewer, labels = lenslab) +
  scale_y_discrete(limits = rev, expand = c(0, 1)) +
  scale_x_continuous(limits = c(0, genomelen), expand = c(0, 1),
                     breaks = chrom_breaks,
                     labels = chrom_labels) +
  theme_bw() +
  labs(x = "ROH position in reference", y = "Individual", color = 'Category') +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 6))

# panel 2: long ROH only
plotzoo2 = zooroh %>% filter(rohcat == '5_Inf')
pp3 <- ggplot(plotzoo2) +
  geom_vline(data = contiglist,
             aes(xintercept = gnPOS),
             color = 'darkgray', linetype = 'dotted') +
  geom_segment(aes(x = gnstart, xend = gnend,
                   y = sample, yend = sample,
                   color = rohcat), size = 3) +
  scale_color_manual(values = rohcatbrewer) +
  scale_y_discrete(limits = rev, expand = c(0, 1)) +
  scale_x_continuous(limits = c(0, genomelen), expand = c(0, 1),
                     breaks = chrom_breaks,
                     labels = chrom_labels) +
  theme_bw() +
  labs(x = "ROH position in reference", y = "Individual", color = 'Category') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 6))

# combine and save --------
pp <- ggarrange(pp2, pp3,
                nrow = 2,
                labels = 'AUTO',
                heights = c(2, 1),
                common.legend = TRUE)

ggsave(filename = 'FigureS8.genomewide_distribution_RZooRoH_newcat_20210912.pdf',
       plot = pp,
       path = plotdir,
       width = 12,
       height = 12)

# write output --------
write.csv(zooroh,
          file = './derived_data/FigS10.csv')

# cleanup --------
closeAllConnections()
