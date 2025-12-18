# Title: Get average ROH lengths and the generation of coalescent (10Mb‐bin grouping)
# Author: Meixi Lin (updated by Kaden Winspear)
# Date: Wed Sep  1 20:36:02 2021 (edited for four‐bin regime)

# ──────────────── Preparation ────────────────
rm(list = ls())
cat("\014")
options(echo = TRUE)

setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/")

library(dplyr)
library(ggplot2)

# ──────────────── Variables ────────────────
today = format(Sys.Date(), "%Y%m%d")

# Define the four ROH length thresholds (in bases) and category labels
rohlens = c(0.1, 1, 5, 10) * 1e+6
rohcats  = c("0.1_1", "1_5", "5_10", "10_Inf")

# ──────────────── Load and recategorize data ────────────────
# Note: this RDS should contain one row per ROH segment with columns 
#       including at least: length (in bases), sample (ID string), and chromosome info.
zooroh = readRDS(file = "./derived_data/ROH_zooroh_summary_final_10mb_update_20250605.rds")

# If “pop” column does not already exist, derive it from “sample” prefix
if (!"pop" %in% colnames(zooroh)) {
  zooroh = zooroh %>%
    mutate(pop = stringr::str_sub(sample, start = 1, end = 3))
}

# Re‐assign each ROH segment into one of four bins:
zooroh = zooroh %>%
  mutate(
    rohcat = case_when(
      length >= rohlens[1] & length <  rohlens[2]  ~ rohcats[1],  # 0.1–1 Mb
      length >= rohlens[2] & length <  rohlens[3]  ~ rohcats[2],  # 1–5 Mb
      length >= rohlens[3] & length <  rohlens[4]  ~ rohcats[3],  # 5–10 Mb
      length >= rohlens[4]                         ~ rohcats[4],  # >10 Mb
      TRUE                                         ~ NA_character_
    )
  ) %>%
  # Drop any segments that failed to fall into a bin (should be none if all lengths ≥ 100 kb)
  filter(!is.na(rohcat)) %>%
  # Convert to factor so that grouping order is consistent
  mutate(rohcat = factor(rohcat, levels = rohcats))

# ──────────────── Compute average lengths per category and population ────────────────
# First, calculate each individual’s mean ROH length within each category:
meanroh = zooroh %>%
  group_by(rohcat, sample, pop) %>%
  summarise(
    meanlen = mean(length),
    .groups = "drop"
  )

# Next, for each population and category, take the average of those individual means:
poproh = meanroh %>%
  group_by(rohcat, pop) %>%
  summarise(
    poplen = mean(meanlen),
    .groups = "drop"
  ) %>%
  # Convert poplen to Mb and compute the inbred‐generation estimate:
  mutate(inbred_gen = 100 / (2 * (poplen / 1e+6)))

# ──────────────── Output ────────────────
write.csv(
  poproh,
  file = paste0("./derived_data/meanroh_zoo_10mb_", today, ".csv"),
  row.names = FALSE
)

# ──────────────── Cleanup ────────────────
date()
closeAllConnections()
