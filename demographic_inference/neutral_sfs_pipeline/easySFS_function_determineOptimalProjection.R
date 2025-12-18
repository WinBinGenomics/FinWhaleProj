# Title: plot the projection preview output 
# Author: Meixi Lin (meixilin@ucla.edu)
# Modified by: Kaden Winspear @ CSUSM to support multiple populations with filenames like ENP_neutral_SNPS_01_PreviewProjection_Rformat.txt
# Date: Updated July 2025

# preparation --------
options(echo = FALSE)
rm(list = ls())
cat("\014")

require(ggplot2)
require(dplyr)

# def variables --------
args = commandArgs(trailingOnly = TRUE)
data.dir = as.character(args[1])
plot.dir = as.character(args[2])
filePREFIX = as.character(args[3])  # should be something like "PreviewProjection"

dir.create(plot.dir, recursive = TRUE)
setwd(data.dir)

# load data --------
fileList <- list.files(
  pattern = paste0(".*_", filePREFIX, "_Rformat\\.txt$"),
  path = data.dir,
  full.names = FALSE
)

# Extract population names from the start of each filename
# Example: ENP_neutral_SNPS_01_PreviewProjection_Rformat.txt → ENP
pops <- sub("_.*$", "", fileList)

# Read and combine all data, tagging with population
alldata <- data.frame()
for (ii in seq_along(fileList)) {
  pop <- pops[ii]
  tempdt <- read.csv(file = fileList[ii]) %>%
    dplyr::mutate(population = pop)
  alldata <- rbind(alldata, tempdt)
}

# main --------
# Get projection with maximum SNPs per population
maxima <- alldata %>%
  group_by(population) %>%
  filter(snps == max(snps))

print("The maxima for projections:")
print(maxima)

# Plot projections
p0 <- ggplot(alldata, aes(x = projection, y = snps)) +
  geom_point() +
  facet_wrap(~population, scales = "free_x") +
  theme_bw() +
  scale_x_continuous(breaks = seq(2, max(alldata$projection), 2))

ggsave(
  filename = paste0(plot.dir, "/", filePREFIX, ".easySFS.projections.allPops.pdf"),
  plot = p0,
  height = 5,
  width = 10
)

# cleanup --------
closeAllConnections()
