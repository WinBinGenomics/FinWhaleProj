# Title: Combine bcftools output
# Author: Meixi Lin
# Date: Wed Jan  6 16:59:51 2021

# preparation --------
rm(list = ls())
cat("\014")

# def functions --------
read_RG_roh <- function(data){
    output <- read.table(file = data, col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"), fill=T, stringsAsFactors = FALSE)
    output1 <- base::subset(output, row_type == "RG")
    return(output1)
}

# def variables --------
outdir = '/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/'
today = format(Sys.Date(), "%Y%m%d")

# load data --------
individuals = read.csv(file = "/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/AdmixRm_fin_popmap.csv", stringsAsFactors = F, header = 1)

# main --------
# read files
rglines = data.frame()
for (ii in 1:nrow(individuals)) {
    sample=individuals[ii,'SampleId']
    pop=individuals[ii,'PopId']
    filename = paste0("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/",
                      pop, "/",sample, "_", pop, "_concat_fwhale_roh_bcftools_G30_ACANGT", ".out.gz")
    if(!file.exists(filename)) {
        print(paste("WARNING:", sample, "No file named", filename))
    } else {
        temp = read_RG_roh(filename)
        rglines = base::rbind(rglines, temp, stringsAsFactors = FALSE)
    }
}

# output lines
write.csv(rglines, file = paste0(outdir, "allindividuals_concat_fwhale_roh_bcftools_G30_ACANGT_", today, ".csv"))

# cleanup --------
closeAllConnections()
