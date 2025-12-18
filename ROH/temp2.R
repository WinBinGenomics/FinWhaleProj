# Date: Wed Sep  8 00:25:56 2021 (with 10Mb-bin updates)
# Related to main text: wilcoxon test and generate data frame for downstream plotting

# ──────────────── Preparation ────────────────
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(ggplot2)
library(reshape2)  # for dcast

setwd('/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/')

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

calculate_nlenrohcat <- function(roh, software) {
  output = roh %>%
    dplyr::group_by(sample, rohcat) %>%
    dplyr::summarise(
      countcat = n(),
      sumcat   = sum(length),
      .groups  = "drop"
    )
  outputn = output %>%
    dcast(sample ~ rohcat, value.var = "countcat")
  colnames(outputn)[-1] = paste0("N_", software, "_", colnames(outputn)[-1])
  outputn[is.na(outputn)] = 0

  outputlen = output %>%
    dcast(sample ~ rohcat, value.var = "sumcat")
  colnames(outputlen)[-1] = paste0("LEN_", software, "_", colnames(outputlen)[-1])
  outputlen[is.na(outputlen)] = 0

  output2 = dplyr::left_join(outputn, outputlen, by = "sample")
  # output both the long and short format
  output1 = output %>%
    dplyr::mutate(software = software)

  return(list(output1, output2))
}

# ──────────────── Variables ────────────────
outdir = "./derived_data/"
dir.create(outdir, showWarnings = FALSE)

today = format(Sys.Date(), "%Y%m%d")

# autosomes total length
genomelen = 2239549461

# ROH length bins (0.1–1, 1–5, 5–10, >10 Mb)
rohlens = c(0.1, 1, 5, 10) * 1e+6
rohcats  = c('0.1_1', '1_5', '5_10', '10_Inf')
lenslab  = c('[0.1, 1) Mb', '[1, 5) Mb', '[5, 10) Mb', '[10, Inf) Mb')

sessionInfo()

# ──────────────── Load Data ────────────────
bcfroh    = readRDS(file = './derived_data/ROH_bcfroh_summary_final_10mb_update_20250605.rds')
zooroh    = readRDS(file = './derived_data/ROH_zooroh_summary_final_10mb_update_20250605.rds')
genomehet = read.csv(
  file           = '/data/shared/snigenda/finwhale_projects/fin_genomics/diversity/genomewide_diversity/autosomalwide_estimates/all70_genomewide_autosomal_heterozygosity_20250404.csv',
  row.names      = 1,
  stringsAsFactors = FALSE
) %>%
  dplyr::select(-starts_with("Total"))

# ──────────────── Main ────────────────

# Get categorized ROH counts and lengths for both methods
bcfroh_nlenrohcat2 = calculate_nlenrohcat(bcfroh, software = "bcf")
zooroh_nlenrohcat2 = calculate_nlenrohcat(zooroh, software = "zoo")

# Long format (for plotting)
forplot = dplyr::bind_rows(
  bcfroh_nlenrohcat2[[1]],
  zooroh_nlenrohcat2[[1]]
)

# Short format (wide counts & lengths)
bcfroh_nlenrohcat = bcfroh_nlenrohcat2[[2]]
zooroh_nlenrohcat = zooroh_nlenrohcat2[[2]]

# Calculate FROH (using only segments ≥ 1 Mb)
bcfroh_froh = calculate_froh(bcfroh, minlen = 1e+6, totallen = genomelen)
zooroh_froh = calculate_froh(zooroh, minlen = 1e+6, totallen = genomelen)

# Merge to create summary table
bcfrohsum = dplyr::left_join(bcfroh_froh, bcfroh_nlenrohcat, by = "sample")
zoorohsum = dplyr::left_join(zooroh_froh, zooroh_nlenrohcat, by = "sample")

rohsum = dplyr::left_join(
  genomehet,
  zoorohsum,
  by = c('SampleId' = 'sample')
) %>%
  dplyr::left_join(
    ., bcfrohsum,
    by = c('SampleId' = 'sample'),
    suffix = c("_zoo", "_bcf")
  )

# Identify samples not called in bcftools
missing_in_bcf = rohsum %>%
  dplyr::filter(is.na(froh_bcf)) %>%
  .$SampleId
# e.g.: "ENPCA09", "ENPOR12", "GOC010"

# Compute totals across all four ROH bins
rohsum2 = rohsum %>%
  dplyr::mutate(
    N_zoo_all  = N_zoo_0.1_1  + N_zoo_1_5  + N_zoo_5_10  + N_zoo_10_Inf,
    N_bcf_all  = N_bcf_0.1_1  + N_bcf_1_5  + N_bcf_5_10  + N_bcf_10_Inf,
    LEN_zoo_all = LEN_zoo_0.1_1 + LEN_zoo_1_5 + LEN_zoo_5_10 + LEN_zoo_10_Inf,
    LEN_bcf_all = LEN_bcf_0.1_1 + LEN_bcf_1_5 + LEN_bcf_5_10 + LEN_bcf_10_Inf
  )

# Run pairwise Wilcoxon test on total zoo counts by population
pairwise.wilcox.test(
  rohsum2$N_zoo_all,
  rohsum2$PopId,
  p.adjust.method = "bonferroni"
)

# ──────────────── Output ────────────────
write.csv(
  rohsum2,
  file = paste0(outdir, "rohsummary_10_Mbzoobcf_", today, ".csv")
)
saveRDS(
  rohsum2,
  file = paste0(outdir, "rohsummary_10MB_zoobcf_", today, ".rds")
)
saveRDS(
  forplot,
  file = paste0(outdir, "rohsummarylong_10Mbzoobcf_", today, ".rds")
)

# ──────────────── Cleanup ────────────────
date()
closeAllConnections()
