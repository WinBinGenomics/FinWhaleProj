# Title: Analyze Fst using the new LDPruned and output tables
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin
# Date: Fri Mar 19 00:42:32 2021

# Modified by Kaden Winspear (CSUSM). For project of analyzing EP fin whale populations. 
# Date -> April 2025

# prep
rm(list = ls())
cat("\014")

library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)

source('/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/plotting_config.R')

# def variables 
today = format(Sys.Date(), "%Y%m%d")

setwd('/Users/kadenwinspear/Documents/finwhale_proj/PopStructure/all70/')
outdir = './derive_data/'
plotdir = './plots/'
logdir = './logs/'

dir.create(outdir)
dir.create(plotdir)
dir.create(logdir)

sink(file = paste0(logdir, 'FST_all70_Blue_', today, '.log'))
dataset = 'all70'
ref = 'Blue'
mafcut = '05'
gdsfile = './data_files/JointCalls_all70_filterpass_biallelic_all_LDPruned_maf05.gds'

sessionInfo()
getwd()

# def functions 
get_allFST <- function(genofile, pop_code) {
    fst_pop <- snpgdsFst(gdsobj = genofile, population=as.factor(pop_code),
                         method="W&C84", autosome.only = F)
    print(paste0('Fst: ', fst_pop$Fst))
    print(paste0('Mean Fst: ',fst_pop$MeanFst))
    return(fst_pop)
}

get_pairFST <- function(genofile, pop.map, popcolname, poporder) {
    pop_code = pop.map[,popcolname]
    fst <- data.frame(t(combn(poporder,2)), stringsAsFactors = FALSE)
    colnames(fst) <- c("Loc1","Loc2")
    fst$Fst <- NA
    fst$MeanFst <- NA
    for(ii in 1:nrow(fst)){
        pair <-  as.character(fst[ii,c("Loc1","Loc2")])
        popid <- which(pop.map[,popcolname] %in% pair)
        sec <- pop.map[popid, c('SampleId', popcolname)]
        colnames(sec) <- c('sample', 'location')
        fstR <- snpgdsFst(gdsobj = genofile, sample.id = sec$sample,population=as.factor(sec$location), method="W&C84",  autosome.only = F, verbose = F, remove.monosnp = T,missing.rate = 0.2)
        fst[ii,3] <- fstR$Fst
        fst[ii,4] <- fstR$MeanFst
    }
    fstc = reshape2::dcast(data = fst, formula = Loc1 ~ Loc2, value.var = 'Fst') %>%
        tibble::column_to_rownames(var = 'Loc1')
    fstc = fstc[poporder[-length(poporder)],poporder[-1]]
    print(round(fstc, digits = 4))
    output = list(fst, fstc)
    return(output)
}

# loads data
genofile = SNPRelate::snpgdsOpen(filename = paste0(ref, '/', gdsfile))

# List of sample ids
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Sample to population identifier. 
popmap <- read.csv(file = "/Users/kadenwinspear/Documents/finwhale_proj/scripts/config/all70_fin_popmap.csv", stringsAsFactors = FALSE)

# Extracts population identifiers
pop_code <- popmap[, 'PopId']

# Compute overall Fst between populations (ESP, ENP, GOC)
fst_pop <- get_allFST(genofile = genofile, pop_code = pop_code)

# Pairwise Fst calculations for ESP, ENP, and GOC
fstpair_pop <- get_pairFST(genofile, popmap, popcolname = 'PopId', poporder = c("ENP", "ESP", "GOC"))

# Output files
write.csv(x = fstpair_pop[[1]], file = paste0(outdir, 'pairFST_pop_maf', mafcut, '_all70_finwhale_', today, '.csv'), row.names = FALSE)
write.csv(x = fstpair_pop[[2]], file = paste0(outdir, 'pairFST_popM_maf', mafcut, '_all70_finwhale_', today, '.csv'), row.names = FALSE)



# cleanup --------
closefn.gds(genofile)
date()
sink()
closeAllConnections()
