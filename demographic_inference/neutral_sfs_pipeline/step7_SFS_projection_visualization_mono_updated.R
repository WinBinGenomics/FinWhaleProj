## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx
## NOTES: R 4.0,v1.0, July 30 2020
# Adjusted by Kaden Winspear @ CSUSM for Eastern Pacific Fin Whale project.
# Plots ESP, ENP, GOC in that order and makes sure each plot has the correct amount of projection bins.
# Fixed monomorphic counts and added proportion plot at the end. 
# July 16 2025

#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Join SFS projections for chromosomes and plot them
#--------------------------------------------------------------------------------------------------

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################

# ------------------------------ Load Required Packages ------------------------------
library(plyr)
library(ggplot2)
library(RColorBrewer)
require(reshape2)

# ------------------------------ Define Functions ------------------------------
Monomorphic <- function(filenames, dataframe, pop){
  c <- 2
  for(file in filenames){
    mono <- read.table(file, header = TRUE)
    dataframe <- cbind(dataframe, mono$HomREFcountPassingProjThreshold[mono$population == pop])
    colnames(dataframe)[c] <- paste0("VCF_", c)
    c <- c + 1
  }
  return(dataframe)
}

Projection <- function(filenames, dataframe, fold){
  c <- 2
  for(file in filenames){
    proj  <- read.table(file, header = TRUE, skip = 1)
    proj2 <- melt(proj)
    dataframe <- cbind(dataframe, proj2[1:fold, 2])  # Keep only 0 to folded bins
    colnames(dataframe)[c] <- paste0("VCF_", c)
    c <- c + 1
  }
  return(dataframe)
}

# ------------------------------ Set Variables ------------------------------
pops <- c("ESP", "ENP", "GOC")
dir <- "/data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral"
dadidir <- file.path(dir, "SFS_projection/SNPs_easySFS_projection_01/dadi")

setwd(dir)

# Monomorphic site files
mono.files <- list.files(paste0(dir, "/SFS_projection_Monomorphic"), pattern = ".perPop", full.names = TRUE)

# ------------------------------ Main Loop ------------------------------
for(pop in pops){
  # --- Monomorphic
  mono <- data.frame(Frequency = 0)
  mono <- Monomorphic(mono.files, mono, pop)
  mono$Counts <- rowSums(mono[, 2:ncol(mono)])
  mono$Pop    <- pop

  # --- SFS projection
  files  <- scan(paste0(dir, "/SFS_projection_fsc/", pop, "_SFS_projection_files.txt"), what = character())
  file1  <- read.table(files[1], skip = 1)
  total_bins <- length(file1)
  folded_bins <- floor(total_bins / 2)

  projection.values <- data.frame(Frequency = 0:folded_bins)
  projection <- Projection(files, projection.values, fold = folded_bins + 1)
  projection$Counts <- rowSums(projection[, 2:ncol(projection)])

  # --- Overwrite monomorphic count correctly
  #Add monomorphic sites
  projection$Counts[1] <- projection$Counts[1] + mono$Counts
  #projection$Counts[1] <- mono$Counts

  # --- Add proportions
  projection$Proportion_all <- projection$Counts / sum(projection$Counts)
  projection$Proportion <- ifelse(projection$Frequency == 0,
                                  NA,
                                  projection$Counts / sum(projection$Counts[projection$Frequency != 0]))
  projection$Pop <- pop

  assign(paste0("projection_", pop), projection)

  # --- Write files
  counts <- projection$Counts
  sfs <- cbind.data.frame(split(counts, rep(1:length(projection$Frequency), times = 1)), stringsAsFactors = FALSE)
  colnames(sfs) <- paste0("d0_", projection$Frequency)

  sink(paste(dir, "/SFS_projection_fsc/", pop, "_MAFpop0.obs", sep = ""))
  cat("1 observation\n")
  write.table(sfs, quote = FALSE, row.names = FALSE)
  sink()

  file <- list.files(dadidir, pattern = paste("^", pop, "-[0-9]+.sfs", sep = ""), full.names = TRUE)
  header <- readLines(file, n = 1)
  sfs_dadi <- read.table(file, skip = 1, header = FALSE)
  sfs_dadi[1, 1:length(projection$Frequency)] <- sfs
  sink(paste(dir, "/SFS_projection_fsc/", pop, "-", as.character(length(sfs_dadi) - 1), ".sfs", sep = ""))
  cat(header, "\n")
  write.table(sfs_dadi, quote = FALSE, row.names = FALSE, col.names = FALSE)
  sink()

  totalSites <- sum(sfs_dadi[1, ])
  sink(paste(dir, "/SFS_projection_fsc/", pop, "-", as.character(length(sfs_dadi) - 1), ".totalSiteCount.L.withMonomorphic.txt", sep = ""))
  cat("pop\ttotalSites\n")
  cat(pop, totalSites, "\n")
  sink()
}

# ------------------------------ Combine & Format for Plotting ------------------------------
projection_allpops1 <- rbind(projection_ESP, projection_ENP, projection_GOC)
projection_allpops2 <- subset(projection_allpops1, Frequency != 0)

projection_allpops1$Pop <- factor(projection_allpops1$Pop, levels = c("ESP", "ENP", "GOC"))
projection_allpops2$Pop <- factor(projection_allpops2$Pop, levels = c("ESP", "ENP", "GOC"))

colors <- c("ESP" = "#7570B3", "ENP" = "#1B9E77", "GOC" = "#D95F02")
axis_text_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

# ------------------------------ Plotting ------------------------------
pdf(paste0(dir, "/SFS_projection_fsc/Final_fwhale_SFS.pdf"))

# 1) Counts, polymorphic only
projection_allpops2$Label <- ave(projection_allpops2$Frequency, projection_allpops2$Pop, FUN = function(x) seq_along(x))
ggplot(projection_allpops2, aes(x = Label, y = Counts, fill = Pop)) +
  geom_bar(stat = "identity") +
  theme_light() + axis_text_theme +
  scale_x_continuous(breaks = projection_allpops2$Label, labels = projection_allpops2$Frequency) +
  labs(title = "Folded SFS projection (poly-only)", x = "Alternative Allele count", y = "Counts") +
  facet_grid(~Pop, scales = "free_x") +
  scale_fill_manual(values = colors) +
  guides(fill = "none")

# 2) Proportions, polymorphic only
projection_allpops2$Label <- ave(projection_allpops2$Frequency, projection_allpops2$Pop, FUN = function(x) seq_along(x))
ggplot(projection_allpops2, aes(x = Label, y = Proportion, fill = Pop)) +
  geom_bar(stat = "identity") +
  theme_light() + axis_text_theme +
  scale_x_continuous(breaks = projection_allpops2$Label, labels = projection_allpops2$Frequency) +
  labs(title = "Folded proportional SFS (poly-only)", x = "Alternative Allele count", y = "Proportion") +
  facet_grid(~Pop, scales = "free_x") +
  scale_fill_manual(values = colors) +
  guides(fill = "none")

# 3) Counts including monomorphic
projection_allpops1$Label <- ave(projection_allpops1$Frequency, projection_allpops1$Pop, FUN = function(x) seq_along(x))
ggplot(projection_allpops1, aes(x = Label, y = Counts, fill = Pop)) +
  geom_bar(stat = "identity") +
  theme_light() + axis_text_theme +
  scale_x_continuous(breaks = projection_allpops1$Label, labels = projection_allpops1$Frequency) +
  labs(title = "Folded SFS projection", subtitle = "Including monomorphic sites",
       x = "Alternative Allele count", y = "Counts") +
  facet_grid(~Pop, scales = "free_x") +
  scale_fill_manual(values = colors) +
  guides(fill = "none")

# 4) Proportions including monomorphic
projection_allpops1$Label <- ave(projection_allpops1$Frequency, projection_allpops1$Pop, FUN = function(x) seq_along(x))
ggplot(projection_allpops1, aes(x = Label, y = Proportion_all, fill = Pop)) +
  geom_bar(stat = "identity") +
  theme_light() + axis_text_theme +
  scale_x_continuous(breaks = projection_allpops1$Label, labels = projection_allpops1$Frequency) +
  labs(title = "Folded proportional SFS", subtitle = "Including monomorphic sites",
       x = "Alternative Allele count", y = "Proportion") +
  facet_grid(~Pop, scales = "free_x") +
  scale_fill_manual(values = colors) +
  guides(fill = "none")

dev.off()
