# Title: Plot admixture
# Author: Meixi Lin
# Adjusted by Kaden Winspear @ CSUSM -> Eastern pacific fin whale project. 

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/plots")
today = format(Sys.Date(), "%Y%m%d")

# sink(file = paste0("logs/admixture_cv_", today,".log"))
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# def functions --------

# def variables --------
mafcut='05'
# mafcut='10'
# mafcut='NA'

sessionInfo()
source("/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/plotting_config.R")

# load data --------
# read cross validation files
adcv = read.csv(file = paste0("/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/Blue/Admixture_20250409/maf", mafcut, "/Admixture_CV_LLsummary_maf", mafcut, ".csv"), stringsAsFactors = FALSE)

forplot = adcv %>%
    dplyr::group_by(K) %>%
    dplyr::summarise(meancv = mean(CVERROR),
                     sdcv = sd(CVERROR))

# main --------
pp1 <- ggplot(forplot, aes(x = K, y = meancv)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=meancv-sdcv, ymax=meancv+sdcv), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = "K", y = "Cross-Validation Error") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 11))

ggsave(filename = paste0("admixture_maf", mafcut, "_cv_", today, ".pdf"), plot = pp1, height = 3, width = 3, bg = "transparent")

write.csv(forplot, file = '/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/derive_data/Admixture_CV.csv')

# cleanup --------
# sink()
closeAllConnections()
