# Title: Mann Whitney U test for the genomewide heterozygosity
# Author: Meixi Li
# Date: Tue Aug 24 20:33:02 2021

# Relevant values: In GOC individuals we found patterns of reduced variation, with an average 1.13 heterozygotes per kb ... In contrast, the ENP population had much higher diversity (1.76 het/kb; Mann-Whitney U [MWU] test p<0.001; Figure 2A)

rm(list = ls())
cat("\014")
options(echo = TRUE)
library(dplyr)
library(reshape2)


# load data --------
totaldt = read.csv(file = '/data/shared/snigenda/finwhale_projects/Summary_stats/all70/Blue/heterozygosity/autosomal/heterozygosityall70_genomewide_autosomal_heterozygosity_20250131.csv', row.names = 1, stringsAsFactors = FALSE)

# main --------
# mwu test
pairwise.wilcox.test(totaldt$GenomeHet, totaldt$PopId)

# Wilcoxon rank sum test
# data:  GenomeHet by PopId
# ENP     ESP    
#ESP 8.6e-11 -      
#GOC 5.2e-13 5.8e-1

# mean values
meandt = totaldt %>%
  dplyr::group_by(PopId) %>%
  dplyr::summarise(PopGenomeHet = mean(GenomeHet),
                   .groups = 'drop')
meandt

# PopId PopGenomeHet
#<chr>        <dbl>
#1 ENP        0.00180
#2 ESP        0.00187
#3 GOC        0.00112
# output files --------

# cleanup --------
date()
closeAllConnections()
