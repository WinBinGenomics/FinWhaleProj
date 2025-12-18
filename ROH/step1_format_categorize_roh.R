# Title: Format and categorize ROH output from both bcftools and zooroh
# Author: Meixi Lin
# Date: Thu Jan  7 09:43:13 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

# require(devtools)
# install_version("RZooRoH", version = "0.2.3", repos = "http://cran.us.r-project.org")
library(RZooRoH)
library(dplyr)
library(ggplot2)

setwd('/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/')

# def functions --------
loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    output = get(ls()[ls() != "fileName"])
    return(output)
}

load_contiglist <- function(contigfile) {
    contiglist = read.csv(file = contigfile, row.names = 1, stringsAsFactors = F)
    # change to zero based and only select the two columns relavant
    contiglist = contiglist %>%
        dplyr::mutate(gnPOS = genomewide_coord -1) %>%
        dplyr::select(SN, gnPOS)
    colnames(contiglist) = c('chrom', 'gnPOS')
    return(contiglist)
}

# format the bcfroh file
format_bcfroh <- function(bcfroh, contiglist) {
    output = dplyr::left_join(bcfroh, contiglist, by = "chrom") %>%
        dplyr::mutate(gnstart = gnPOS + start,
                      gnend = gnPOS + end)
    return(output)
}


# format zooroh
# chrom marked from 1 to 21
format_zooroh <- function(zres, contiglist) {
    sampleids = data.frame(unique(zres@hbdseg$id), zres@sampleids)
    colnames(sampleids) = c("id", "sample")
    contiglist = contiglist %>%
        dplyr::select(-chrom) %>%
        tibble::rowid_to_column(var = "chrom")
    hbdseg = zres@hbdseg
    hbdseg = dplyr::left_join(hbdseg, y = sampleids, by = "id") %>%
        dplyr::left_join(., y = contiglist, by = "chrom") %>%
        dplyr::mutate(gnstart = gnPOS + start_pos,
                      gnend = gnPOS + end_pos) %>%
        dplyr::select(-id)
    return(hbdseg)
}

# categorize roh into three categories
categorize_roh <- function(rohdt, rohlens, rohcats) {
    output = rohdt %>%
        dplyr::mutate(rohcat = case_when(
            length >= rohlens[1] & length < rohlens[2] ~ rohcats[1],
            length >= rohlens[2] & length < rohlens[3] ~ rohcats[2],
            length >= rohlens[3] ~ rohcats[3],
            TRUE ~ "FAIL")) %>%
        dplyr::filter(rohcat != "FAIL") %>%
        dplyr::mutate(pop = stringr::str_sub(sample, start = 1, end = 3))
    return(output)
}

# def variables --------
outdir = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/"
dir.create(outdir)

today = format(Sys.Date(), "%Y%m%d")
#genomelen = 2324429847

#autosomes length
genomelen = 2239549461 

# new roh length and categories (August 2021)
rohlens = c(0.1,1,5)*1e+6
rohcats = c('0.1_1','1_5','5_Inf')
lenslab <- c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb')

sessionInfo()
# source("./scripts/config/plotting_config.R")

# load data --------
bcfroh = read.csv(file="/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/derived_data/allindividuals_concat_fwhale_roh_bcftools_G30_ACANGT_20250602.csv", stringsAsFactors = F, row.names = 1)
zoorohenp = loadRData(fileName = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH/mix10R_base3/ENP_ZooROH_Results.RData")
zoorohgoc = loadRData(fileName = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH/mix10R_base3/GOC_ZooROH_Results.RData")
zoorohesp = loadRData(fileName = "/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/RZOOROH/mix10R_base3/ESP_ZooROH_Results.RData")


contiglist = load_contiglist(contigfile = "/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/blue_contig_summary_with_coord.csv")

# main --------
# format bcfroh ========
bcfroh = format_bcfroh(bcfroh, contiglist)
bcfroh = categorize_roh(bcfroh, rohlens, rohcats)

# format zooroh ========

zooroh = base::rbind(format_zooroh(zres = zoorohenp, contiglist),
                     format_zooroh(zres = zoorohgoc, contiglist),
                     format_zooroh(zres = zoorohesp, contiglist))


zooroh = categorize_roh(zooroh, rohlens, rohcats)

# write data ========
saveRDS(object = zooroh, file = paste0(outdir, "ROH_zooroh_summary_final_", today, ".rds"))
saveRDS(object = bcfroh, file = paste0(outdir, "ROH_bcfroh_summary_final_", today, ".rds"))

# cleanup --------
closeAllConnections()
