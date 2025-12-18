# Title: Format and categorize ROH output from both bcftools and zooroh (10Mb update)
# Author: Meixi Lin (edited by Kaden Winspear)
# Date: Thu Jan  7 09:43:13 2021 (with 10Mb bin updates)
#
# This script now categorizes ROH segments into four bins:
#   0.1–1 Mb, 1–5 Mb, 5–10 Mb, and >10 Mb.
# The resulting RDS filenames include “10mb_update” to distinguish this version.

# ──────────────── Preparation ────────────────
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(RZooRoH)
library(dplyr)
library(ggplot2)
# (stringr::str_sub is used explicitly, so no need to library(stringr) here)

setwd('/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/')

# ──────────────── Helper Functions ────────────────

loadRData <- function(fileName){
  # loads an RData file, and returns the main object inside it
  load(fileName)
  output = get(ls()[ls() != "fileName"])
  return(output)
}

load_contiglist <- function(contigfile) {
  contiglist = read.csv(file = contigfile, row.names = 1, stringsAsFactors = FALSE)
  # change to zero-based coordinate and select only the two relevant columns
  contiglist = contiglist %>%
    dplyr::mutate(gnPOS = genomewide_coord - 1) %>%
    dplyr::select(SN, gnPOS)
  colnames(contiglist) = c('chrom', 'gnPOS')
  return(contiglist)
}

format_bcfroh <- function(bcfroh, contiglist) {
  output = dplyr::left_join(bcfroh, contiglist, by = "chrom") %>%
    dplyr::mutate(
      gnstart = gnPOS + start,
      gnend   = gnPOS + end
    )
  return(output)
}

format_zooroh <- function(zres, contiglist) {
  # Convert ZooROH output (zres) into a data.frame with genome-wide start/end coords
  sampleids = data.frame(unique(zres@hbdseg$id), zres@sampleids)
  colnames(sampleids) = c("id", "sample")

  contiglist = contiglist %>%
    dplyr::select(-chrom) %>%
    tibble::rowid_to_column(var = "chrom")

  hbdseg = zres@hbdseg
  hbdseg = dplyr::left_join(hbdseg, y = sampleids, by = "id") %>%
    dplyr::left_join(., y = contiglist, by = "chrom") %>%
    dplyr::mutate(
      gnstart = gnPOS + start_pos,
      gnend   = gnPOS + end_pos
    ) %>%
    dplyr::select(-id)

  return(hbdseg)
}

# ──────────────── Categorization Function (FOUR bins) ────────────────

categorize_roh <- function(rohdt, rohlens, rohcats) {
  output = rohdt %>%
    dplyr::mutate(rohcat = case_when(
      length >= rohlens[1] & length <  rohlens[2]  ~ rohcats[1],  # 0.1–1 Mb
      length >= rohlens[2] & length <  rohlens[3]  ~ rohcats[2],  # 1–5 Mb
      length >= rohlens[3] & length <  rohlens[4]  ~ rohcats[3],  # 5–10 Mb
      length >= rohlens[4]                         ~ rohcats[4],  # >10 Mb
      TRUE                                         ~ "FAIL"
    )) %>%
    dplyr::filter(rohcat != "FAIL") %>%
    dplyr::mutate(pop = stringr::str_sub(sample, start = 1, end = 3))
  return(output)
}

# ──────────────── Main Variables ────────────────

outdir = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/"
dir.create(outdir, showWarnings = FALSE)

today = format(Sys.Date(), "%Y%m%d")

# autosomes total length (not directly used below, but kept for reference)
genomelen = 2239549461

# Define four ROH-length bins (0.1, 1, 5, 10 Mb) and category labels
rohlens = c(0.1, 1, 5, 10) * 1e+6
rohcats  = c('0.1_1', '1_5', '5_10', '10_Inf')
lenslab  = c('[0.1, 1) Mb', '[1, 5) Mb', '[5, 10) Mb', '[10, Inf) Mb')

sessionInfo()

# ──────────────── Load Input Data ────────────────

# BCFTools ROH summary (has columns “chrom”, “start”, “end”, “length”, “sample”, etc.)
bcfroh = read.csv(
  file      = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/derived_data/allindividuals_concat_fwhale_roh_bcftools_G30_ACANGT_20250602.csv",
  stringsAsFactors = FALSE,
  row.names = 1
)

# ZooROH results for each population
zoorohenp = loadRData(fileName = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH/mix10R_base3/ENP_ZooROH_Results.RData")
zoorohgoc = loadRData(fileName = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH/mix10R_base3/GOC_ZooROH_Results.RData")
zoorohesp = loadRData(fileName = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH/mix10R_base3/ESP_ZooROH_Results.RData")

# Contig list for translating chromosome IDs into genome-wide positions
contiglist = load_contiglist(
  contigfile = "/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/blue_contig_summary_with_coord.csv"
)

# ──────────────── Format & Categorize BCFTools ROH ────────────────

bcfroh = format_bcfroh(bcfroh, contiglist)
bcfroh = categorize_roh(bcfroh, rohlens, rohcats)

# ──────────────── Format & Categorize ZooROH Results ────────────────

zooroh = base::rbind(
  format_zooroh(zres = zoorohenp, contiglist),
  format_zooroh(zres = zoorohgoc, contiglist),
  format_zooroh(zres = zoorohesp, contiglist)
)

zooroh = categorize_roh(zooroh, rohlens, rohcats)

# ──────────────── Save Outputs (with “10mb_update” in filenames) ────────────────

saveRDS(
  object = zooroh,
  file   = paste0(outdir, "ROH_zooroh_summary_final_10mb_update_", today, ".rds")
)

saveRDS(
  object = bcfroh,
  file   = paste0(outdir, "ROH_bcfroh_summary_final_10mb_update_", today, ".rds")
)

# ──────────────── Cleanup ────────────────
closeAllConnections()
