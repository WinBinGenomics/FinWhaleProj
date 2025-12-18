# Modification: Change to the final directory (with 10Mb-bin updates)
# Date: Sun Sep 12 10:43:27 2021 (edited by Kaden Winspear)

# ──────────────── Preparation ────────────────
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/")

library(dplyr)
library(ggplot2)

# ──────────────── Functions ────────────────

calculate_froh <- function(roh, minlen = 1e+6, totallen) {
  output = roh %>%
    dplyr::filter(length >= minlen) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      froh = sum(length) / totallen,
      .groups = "drop"
    )
  return(output)
}

join_frohdtl <- function(frohdtl, frohnames) {
  # Assumes frohdtl is a list of data.frames, each with columns (sample, froh)
  frohdt = dplyr::full_join(x = frohdtl[[1]], y = frohdtl[[2]], by = "sample")
  frohdt = dplyr::full_join(x = frohdt,    y = frohdtl[[3]], by = "sample")
  frohdt = dplyr::full_join(x = frohdt,    y = frohdtl[[4]], by = "sample")
  colnames(frohdt) = frohnames
  return(frohdt)
}

# ──────────────── Variables ────────────────
indir  = "./derived_data/"
outdir = "./derived_data/"

today = format(Sys.Date(), "%Y%m%d")

# autosomes total length
genomelen = 2239549461

# ROH length thresholds (0.1 Mb, 1 Mb, 5 Mb, 10 Mb)
rohlens = c(0.1, 1, 5, 10) * 1e+6
rohcats = c('0.1_1', '1_5', '5_10', '10_Inf')  # (not directly used here)
lenslab = c('[0.1, 1) Mb', '[1, 5) Mb', '[5, 10) Mb', '[10, Inf) Mb')

sessionInfo()

# ──────────────── Load Data ────────────────
bcfroh = readRDS(file = './derived_data/ROH_bcfroh_summary_final_10mb_update_20250605.rds')
zooroh = readRDS(file = './derived_data/ROH_zooroh_summary_final_10mb_update_20250605.rds')

# ──────────────── Main ────────────────

# Calculate FROH at each threshold for ZooROH
# (List of four data.frames: each contains sample and froh for >= rohlens[i])
frohdtl_zoo = lapply(
  X = rohlens,
  FUN = function(th) calculate_froh(roh = zooroh, minlen = th, totallen = genomelen)
)

# Calculate FROH at each threshold for BCFTools ROH
frohdtl_bcf = lapply(
  X = rohlens,
  FUN = function(th) calculate_froh(roh = bcfroh, minlen = th, totallen = genomelen)
)

# Join into a single data.frame with four columns + sample
# Names correspond to thresholds: 0.1Mb (100k), 1Mb, 5Mb, 10Mb
frohdt_zoo = join_frohdtl(
  frohdtl   = frohdtl_zoo,
  frohnames = c(
    "Sample",
    "F_ROH_100k_zoo",
    "F_ROH_1M_zoo",
    "F_ROH_5M_zoo",
    "F_ROH_10M_zoo"
  )
)

frohdt_bcf = join_frohdtl(
  frohdtl   = frohdtl_bcf,
  frohnames = c(
    "Sample",
    "F_ROH_100k_bcf",
    "F_ROH_1M_bcf",
    "F_ROH_5M_bcf",
    "F_ROH_10M_bcf"
  )
)

# Merge ZooROH and BCFTools results into one table
froh_all = dplyr::full_join(frohdt_zoo, frohdt_bcf, by = "Sample")

# ──────────────── Output Files ────────────────
write.csv(
  froh_all,
  file = paste0(outdir, "FROH_100k1M5M10Mb_all_", today, ".csv"),
  row.names = FALSE,
  quote = FALSE
)

# ──────────────── Cleanup ────────────────
date()
closeAllConnections()
