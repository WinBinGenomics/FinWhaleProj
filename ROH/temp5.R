options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH")

library(dplyr)
library(ggplot2)
library(ggpubr)

# mapping scaffold IDs to chromosome numbers --------
chromosome_map <- c(
  "NC_045785.1" = "1",  "NC_045786.1" = "2",  "NC_045787.1" = "3",
  "NC_045788.1" = "4",  "NC_045789.1" = "5",  "NC_045790.1" = "6",
  "NC_045791.1" = "7",  "NC_045792.1" = "8",  "NC_045793.1" = "9",
  "NC_045794.1" = "10", "NC_045795.1" = "11", "NC_045796.1" = "12",
  "NC_045797.1" = "13", "NC_045798.1" = "14", "NC_045799.1" = "15",
  "NC_045800.1" = "16", "NC_045801.1" = "17", "NC_045802.1" = "18",
  "NC_045803.1" = "19", "NC_045804.1" = "20", "NC_045805.1" = "21"
)

# def functions --------
load_contiglist <- function(contigfile) {
  contiglist <- read.csv(file = contigfile, row.names = 1, stringsAsFactors = FALSE)
  contiglist <- contiglist %>%
    dplyr::mutate(
      gnPOS = genomewide_coord - 1,
      chrom = chromosome_map[as.character(SN)]
    ) %>%
    dplyr::select(chrom, gnPOS)
  colnames(contiglist) <- c("chrom", "gnPOS")
  return(contiglist)
}

# def variables --------
indir    <- "./derived_data/"
outdir   <- "./derived_data/"
plotdir  <- "./plots/"
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

today     <- format(Sys.Date(), "%Y%m%d")
genomelen <- 2239549461  # total autosomal length

# load data --------
zooroh <- readRDS(file = "./derived_data/ROH_zooroh_summary_final_10mb_update_20250605.rds")
contiglist <- load_contiglist(
  contigfile = "/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/blue_contig_summary_with_coord.csv"
)

# ensure ROH categories are factors in the correct order --------
rohlens <- c(0.1, 1, 5, 10) * 1e+6
rohcats <- c("0.1_1", "1_5", "5_10", "10_Inf")
zooroh <- zooroh %>%
  mutate(
    rohcat = factor(rohcat, levels = rohcats)
  )

# 1) Extract and sort sample IDs by population prefix --------
esp_ids <- sort(unique(zooroh$sample[grepl("^ESP", zooroh$sample)]))
enp_ids <- sort(unique(zooroh$sample[grepl("^ENP", zooroh$sample)]))
goc_ids <- sort(unique(zooroh$sample[grepl("^GOC", zooroh$sample)]))

# 2) Concatenate in the original desired order: ESP → ENP → GOC
orig_order_ids <- c(esp_ids, enp_ids, goc_ids)

# 3) Reverse so ESP ends up on top
rev_ids <- rev(orig_order_ids)

# 4) Convert 'sample' column into a factor with these reversed levels
zooroh <- zooroh %>%
  mutate(sample = factor(sample, levels = rev_ids))

# prepare axis breaks and labels --------
chrom_breaks <- contiglist$gnPOS
chrom_labels <- contiglist$chrom

# readable labels for legend --------
lenslab <- c("[0.1, 1) Mb", "[1, 5) Mb", "[5, 10) Mb", "[10, Inf) Mb")

sessionInfo()

# main plotting --------
library(paletteer)

# Get the last 4 colors from the Blue2DarkOrange18Steps palette
cb_colors <- paletteer_d("colorBlindness::Blue2DarkOrange18Steps")[15:18]

col_0to1   <- cb_colors[1]  # for 0.1–1 Mb
col_1to5   <- cb_colors[2]  # for 1–5 Mb
col_5to10  <- cb_colors[3]  # for 5–10 Mb
col_10inf  <- cb_colors[4]  # for >10 Mb



# Combine into a named vector in the order of rohcats
cols <- c(
  col_0to1,    # 0.1–1 Mb
  col_1to5,    # 1–5 Mb
  col_5to10,   # 5–10 Mb
  col_10inf    # >10 Mb
)
names(cols) <- rohcats

# panel 1: all ROH segments
pp2 <- ggplot(zooroh) +
  geom_segment(
    aes(
      x = gnstart, xend = gnend,
      y = sample, yend = sample,
      color = rohcat
    ),
    size = 3
  ) +
  scale_color_manual(values = cols, labels = lenslab) +
  scale_y_discrete(expand = c(0, 1)) +
  scale_x_continuous(
    limits = c(0, genomelen),
    expand = c(0, 1),
    breaks = chrom_breaks,
    labels = chrom_labels
  ) +
  theme_bw() +
  labs(
    x = "Chromosomes",
    y = "Individual",
    color = "Category"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 0,
      vjust = 1,
      hjust = 0.5,
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 6)
  )

# panel 2: ROH categories 5–10 Mb AND >10 Mb
plotzoo2 <- zooroh %>% filter(rohcat %in% c("5_10", "10_Inf"))

pp3 <- ggplot(plotzoo2) +
  geom_vline(
    data = contiglist,
    aes(xintercept = gnPOS),
    color = "darkgray",
    linetype = "dotted"
  ) +
  geom_segment(
    aes(
      x = gnstart, xend = gnend,
      y = sample, yend = sample,
      color = rohcat
    ),
    size = 3
  ) +
  scale_color_manual(values = cols, labels = lenslab) +
  scale_y_discrete(expand = c(0, 1)) +
  scale_x_continuous(
    limits = c(0, genomelen),
    expand = c(0, 1),
    breaks = chrom_breaks,
    labels = chrom_labels
  ) +
  theme_bw() +
  labs(
    x = "Chromosomes",
    y = "Individual",
    color = "Category"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      vjust = 1,
      hjust = 0.5,
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 6)
  )

# combine and save --------
pp <- ggarrange(
  pp2, pp3,
  nrow = 2,
  labels = "AUTO",
  heights = c(2, 1),
  common.legend = TRUE
)

ggsave(
  filename = "genomewide_distribution_RZooRoH_newcat_20210912.png",
  plot = pp,
  path = plotdir,
  width = 12,
  height = 12
)

# write output --------
write.csv(
  zooroh,
  file = "./derived_data/FigS10Mb.csv",
  row.names = FALSE,
  quote = FALSE
)

# cleanup --------
closeAllConnections()
