##### Moreno Lab
## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx
## NOTES: R 4.0,v1.0, July 23 2020
# adjusted by Kaden Winspear @ CSUSM - > For eastern pacific fin whales. Added small additions to allow three pops to be plotted + color.

#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Put together the preview of SFS projections from each chromosomes/scaffold
#--------------------------------------------------------------------------------------------------

# Usage in hoffman: Rscript SFS_preview_v1.R

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################

# Load the R packages
library(plyr)
library(readr)
library(ggplot2)
library(RColorBrewer)

#Define functions ------------------------------
readtables <- function(filenames, dataframe){
  c <- 1
  for(file in filenames){
    dataframe <- cbind(dataframe, read_csv(file, skip = 1, col_names = c("prj", paste0("SNPs_", c)), show_col_types = FALSE)[2])
    c <- c+1
  }
  return(dataframe)
}

#Define variables --------------------
dir <- "/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SNPs"
dir.create("/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/Preview", showWarnings = FALSE)
outputdir <- "/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/Preview"
pops <- c("ESP", "ENP", "GOC")

#Main---------------------

#Obtain tables for each pop that contain projection values of each chromosome file
for(pop in pops){
  preview.files <- list.files(dir, pattern = pop, full.names = TRUE)
  preview.values <- read_csv(preview.files[1], col_names = c("Projection","snp"), skip = 1, show_col_types = FALSE)[1]
  preview <- readtables(preview.files, preview.values)
  write.table(preview, paste0(outputdir,"/SFS_preview_",pop,"_Sergio.txt"), quote = FALSE, row.names = FALSE)
  preview$SNPs <- rowSums(preview[ , -1])
  preview$pop <- pop
  assign(paste0("preview_", pop), preview)
}

#Join populations
preview_allpops <- rbind(preview_ESP, preview_ENP, preview_GOC)

# Custom colors: ESP = purple, ENP = green, GOC = orange
custom_colors <- c("ENP" = "#1b9e77", "ESP" = "#984ea3", "GOC" = "#ff7f00")

# Set pop order for faceting
preview_allpops$pop <- factor(preview_allpops$pop, levels = c("ESP", "ENP", "GOC"))

# Plot preview
p1 <- ggplot(preview_allpops, aes(x = Projection, y = SNPs, group = pop, color = pop)) +
  geom_line() + geom_point() +
  facet_wrap(~pop, scales = "free_x") +
  theme_light() +
  scale_x_continuous(breaks = seq(2, 80, 2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5)) +
  labs(title = "Site Frequency Spectrum Projection Preview", y = "Number of segregating sites") +
  guides(color = FALSE) +
  scale_color_manual(values = custom_colors)


ggsave(paste0(outputdir, "/FinWhale_easySFS_projections_allPops.pdf"), p1, height = 5, width = 10)
