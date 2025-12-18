## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx
## NOTES: R 4.0,v1.0, July 30 2020
# adjusted by Kaden Winspear @ CSUSM for Eastern Pacific Fin Whale project.
# Plots ESP, ENP, GOC in that order and makes sure each plot has the correct amount of projection bins.
# July 16 2025

#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Join SFS projections for chromosomes and plot them
#--------------------------------------------------------------------------------------------------

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################



# Load packages
library(ggplot2)
library(RColorBrewer)
require(reshape2)

# Define function to read projection tables
readtables <- function(filenames, dataframe, fold) {
  c <- 1
  for(file in filenames){
    proj <- read.table(file, header = TRUE, skip = 1)
    proj2 <- melt(proj)
    dataframe <- cbind(dataframe, proj2[2:fold, 2])
    colnames(dataframe)[c+1] <- paste0("VCF_", c)
    c <- c + 1
  }
  return(dataframe)
}

# Set directory and populations
dir <- "/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection_fsc"
pops <- c("ESP", "ENP", "GOC")  # Desired plot order

# Set working directory
setwd(dir)

# Loop through each population
for(pop in pops){
  files <- scan(paste0(pop, "_SFS_projection_files.txt"), what = character())
  file1 <- read.table(files[1], skip = 1)
  ssfold <- (length(file1) + 1)/2
  projection.values <- data.frame("Frequency" = 2:ssfold - 1)
  projection <- readtables(files, projection.values, ssfold)
  projection$Counts <- rowSums(projection[, 2:ncol(projection)])
  write.table(projection, paste0("SFS_projection_", pop), quote = FALSE, row.names = FALSE)
  projection$Proportion <- projection$Counts / sum(projection$Counts)
  projection$Pop <- pop
  assign(paste0("projection_", pop), projection)
}

# Combine populations
projection_allpops <- rbind(projection_ESP, projection_ENP, projection_GOC)
projection_allpops$Pop <- factor(projection_allpops$Pop, levels = pops)

# Manually assign Dark2 colors in desired order
colors <- c("ESP" = "#7570B3",  # purple
            "ENP" = "#1B9E77",  # green
            "GOC" = "#D95F02")  # orange

# Plot
pdf("SFS.pdf")

# Absolute SFS
ggplot(projection_allpops, aes(x = Frequency, y = Counts, fill = Pop)) +
  geom_bar(stat = "identity") +
  theme_light() +
  scale_x_continuous(breaks = unique(projection_allpops$Frequency)) +
  labs(title = "Folded SFS projection", x = "Alternative Allele count", y = "Counts") +
  facet_grid(~Pop, scales = "free_x") +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

# Proportional SFS
ggplot(projection_allpops, aes(x = Frequency, y = Proportion, fill = Pop)) +
  geom_bar(stat = "identity") +
  theme_light() +
  scale_x_continuous(breaks = unique(projection_allpops$Frequency)) +
  labs(title = "Folded proportional SFS", x = "Alternative Allele count", y = "Proportion") +
  facet_grid(~Pop, scales = "free_x") +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

dev.off()
