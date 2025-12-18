# Title: configuration for plot
# Author: kaden winspear

# def variables --------
# plotting colors
library(RColorBrewer)

mypops = c("ENP", "GOC", "ESP")

dark2_colors = brewer.pal(n = 3, name = "Dark2")

mycolors = c("ENP" = dark2_colors[1],  # Green
             "GOC" = dark2_colors[2],  # Orange
             "ESP" = dark2_colors[3])  # Purple

names(mycolors) = mypops

mylocs = c("ENP","GOC","ESP")
loccolors = c("ENP" = dark2_colors[1],  # Green
             "GOC" = dark2_colors[2],  # Orange
             "ESP" = dark2_colors[3])  # Purple

names(loccolors) = mylocs
poporder = c("ENP", "GOC", "ESP")

myspecies = c("BalAcu", "BalMus", "EubGla", "MegNov", "BalPhy")
speccolors = c("#F0A0FF","#0075DC","#4C005C","#94FFB5","#808080")

names(speccolors) = myspecies

# pals::pal.bands(speccolors)

# discarding individuals
# ENPOR12: low genotype depth
# GOC077: low genotype depth 
# ENPCA01: Admixture individuals
# ENPCA09: Admixture individuals
# GOC010: Admixture individuals
# GOC080: Related to GOC006 in relateness analyses
# GOC111: Related to GOC025 in relateness analyses
# GOC053: Related to GOC071 in relateness analyses
discardind = c("ENPOR12", "ENPCA01", "ENPCA09", "GOC010", "GOC080", "GOC111", "GOC077","GOC053")
