# Title: Output summaries for bcftools and zooroh comparisons
# Author: Meixi Lin

# Adjusted by Kaden Winspear to Include Eastern South Pacific Population.
# June 2025

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/")

library(dplyr)
library(ggplot2)
library(ggpmisc)
library(reshape2)
library(tidyr)

# def functions --------
get_rohnames <- function(rohlens) {
    output = character(0)
    for (ii in 1:length(rohlens)) {
        if (ii == length(rohlens)) {
            temp = paste0(rohlens[ii], '_Inf')
        } else {
            temp = paste0(rohlens[ii], '_', rohlens[ii + 1])
        }
        output = c(output, temp)
    }
    return(output)
}

get_forplot <- function(rohsum) {
    forplot = rohsum %>%
        select(-starts_with('froh'), -GenomeHet) %>%
        reshape2::melt(., id.vars = c('SampleId', 'PopId'))
    dt = reshape2::colsplit(forplot$variable, pattern = '_', names = c('type', 'software', 'length'))
    forplot = cbind(forplot, dt)
    forplotbcf = forplot %>%
        filter(software == 'bcf') %>%
        rename(bcftools = value) %>%
        select(SampleId, PopId, type, length, bcftools)
    forplotzoo = forplot %>%
        filter(software == 'zoo') %>%
        rename(RZooRoH = value) %>%
        select(SampleId, PopId, type, length, RZooRoH)
    forplot = full_join(x = forplotbcf, y = forplotzoo,
                        by = c('SampleId', 'PopId', 'type', 'length')) %>%
        drop_na() %>%
        mutate(type = ifelse(type == 'N', 'Total number', 'Total length (Mb)'))
    return(forplot)
}

# def variables --------
plotdir = "./plots/"
outdir = "./derived_data/"
dir.create(plotdir, recursive = TRUE)
dir.create(outdir, recursive = TRUE)

today = format(Sys.Date(), "%Y%m%d")
genomelen = 2239549461

rohlens = c(0.1, 1, 5) * 1e+6
rohcats = get_rohnames(c(0.1, 1, 5))

sessionInfo()
source("/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/plotting_config.R")

# load data --------
rohsum2 = readRDS(file = './derived_data/rohsummary_zoobcf_20250603.rds')

# plot comparisons ========
forplot <- get_forplot(rohsum2)
types <- unique(forplot$type)
lens <- unique(forplot$length)

lenslab <- c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb', '[0.1, Inf) Mb')
lenscat <- c('short', 'medium', 'long', 'all')
plotlist <- vector(mode = "list", length = length(types) * length(lens))
counter = 0

for (ii in types) {
    for (j in 1:length(lens)) {
        jj = lens[j]
        lenlab = lenslab[j]
        lencat = lenscat[j]
        counter = counter + 1
        forplotx = forplot %>%
            filter(type == ii, length == jj)
        if (ii == "Total length (Mb)") {
            forplotx = forplotx %>%
                mutate(bcftools = bcftools / 1e+6,
                       RZooRoH = RZooRoH / 1e+6)
        }
        pltrange = range(forplotx$bcftools, forplotx$RZooRoH)
        my.formula = y ~ x

        pp <- ggplot() +
            # raw points (bcftools vs. RZooRoH)
            geom_point(
                data = forplotx,
                aes(x = bcftools, y = RZooRoH),
                shape = "."
            ) +
            # overall (black dashed) linear fit
            geom_smooth(
                data = forplotx,
                aes(x = bcftools, y = RZooRoH),
                method = "lm",
                se = FALSE,
                formula = my.formula,
                size = 0.5,
                linetype = "dashed",
                color = "black"
            ) +
            # equation + R² for the overall fit (default placement)
            ggpmisc::stat_poly_eq(
                data = forplotx,
                mapping = aes(
                    x = bcftools,
                    y = RZooRoH,
                    label = paste(..eq.label.., ..rr.label.., sep = "~~~")
                ),
                formula = my.formula,
                parse = TRUE,
                size = 3
            ) +
            # per‐population colored linear fits
            geom_smooth(
                data = forplotx,
                aes(x = bcftools, y = RZooRoH, color = PopId),
                method = "lm",
                se = FALSE,
                formula = my.formula,
                size = 0.8
            ) +
            # raw colored points by PopId
            geom_point(
                data = forplotx,
                aes(x = bcftools, y = RZooRoH, color = PopId)
            ) +
            # per‐population equations, lowered to npc.y = 0.60
            ggpmisc::stat_poly_eq(
                data = forplotx,
                mapping = aes(
                    x = bcftools,
                    y = RZooRoH,
                    color = PopId,
                    label = paste(..eq.label.., ..rr.label.., sep = "~~~")
                ),
                formula = my.formula,
                parse = TRUE,
                size = 3,
                geom = "text_npc",
                npc.x = "right",   # keep them right‐aligned
                npc.y = 0.60,      # moved down to 60% height
                vstep = 0.08
            ) +
            # 1:1 reference line
            geom_abline(
                slope = 1,
                intercept = 0,
                color = "darkgray",
                linetype = "dotted"
            ) +
            labs(
                title = paste(ii, "of", lencat, "ROH"),
                subtitle = paste("ROH length:", lenlab)
            ) +
            coord_fixed(
                ratio = 1,
                xlim = pltrange,
                ylim = pltrange
            ) +
            scale_color_manual(values = mycolors) +
            theme_bw() +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 12),
                plot.subtitle = element_text(size = 10)
            )

        plotlist[[counter]] <- pp
    }
}

ppout <- ggpubr::ggarrange(
    plotlist = plotlist,
    nrow = length(types),
    ncol = length(lens)
)

# output --------
ggsave(
    filename = paste0("FigureS8.ROH_bcfzoo_compare_", today, ".pdf"),
    path = plotdir,
    width = 14,
    height = 8
)

write.csv(forplot, file = "./derived_data/FigS9.csv")

# cleanup --------
date()
closeAllConnections()
