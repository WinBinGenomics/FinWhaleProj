# Dr. Robinson UCLA.
# This script was first written by Dr. Robinson and was adjusted / modified for this project.
# Adjusted / modified -> Kaden Winspear @ CSUSM Eastern Pacific Fin Whale Project 
# Date April 2025

library(plyr)  # Plyr | splits, manipulates, and recombines data.

# ──────────────── Setup ────────────────
setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/diversity/sliding_window/1mbwins")
all_files      <- list.files(pattern = "^combined_.*\\.txt$", ignore.case = TRUE)
search_strings <- c("ENPAK19","ENPAK20","ENPAK21","ENPAK22","ENPAK23","ENPAK24","ENPAK25","ENPAK26","ENPAK27","ENPAK28","ENPAK29","ENPAK30","ENPBC16","ENPBC17","ENPBC18","ENPCA01","ENPCA02","ENPCA03","ENPCA04","ENPCA05","ENPCA06","ENPCA07","ENPCA08","ENPCA09",
                    "ENPOR10","ENPOR11","ENPOR12","ENPOR13","ENPWA14","ENPWA15","ESPCL01","ESPCL02","ESPCL03","ESPCL04","ESPCL05","ESPCL06","ESPCL07","ESPCL08","ESPCL09","ESPCL10","ESPCL11","ESPCL12","ESPCL13","ESPCL14","ESPCL15","ESPCL16","ESPCL17","ESPCL18","ESPCL19","ESPCL20",
                    "GOC002","GOC006","GOC010","GOC025","GOC038","GOC050","GOC053","GOC063","GOC068","GOC071","GOC077","GOC080","GOC082","GOC086","GOC091","GOC100","GOC111","GOC112","GOC116","GOC125")
winsize <- 1e6

# ──────────────── Functions ────────────────
plotwinhet <- function(samplename, colset, temp) {
  calls_col <- paste0("calls_", samplename)
  hets_col  <- paste0("hets_",  samplename)
  temp[[calls_col]] <- as.numeric(temp[[calls_col]])
  temp[[hets_col]]  <- as.numeric(temp[[hets_col]])
  meanhet <- sum(temp[[hets_col]], na.rm = TRUE) / sum(temp[[calls_col]], na.rm = TRUE)
  plotname <- paste(samplename, "\nmean het. = ", sprintf("%.3f", 1000 * meanhet), " per kb", sep = "")
  
  chromosome_map <- c(
    "NC_045785.1"="1","NC_045786.1"="2","NC_045787.1"="3","NC_045788.1"="4","NC_045789.1"="5",
    "NC_045790.1"="6","NC_045791.1"="7","NC_045792.1"="8","NC_045793.1"="9","NC_045794.1"="10",
    "NC_045795.1"="11","NC_045796.1"="12","NC_045797.1"="13","NC_045798.1"="14","NC_045799.1"="15",
    "NC_045800.1"="16","NC_045801.1"="17","NC_045802.1"="18","NC_045803.1"="19","NC_045804.1"="20",
    "NC_045805.1"="21"
  )
  temp$chrom <- mapvalues(temp$chrom, from = names(chromosome_map), to = chromosome_map)
  
  pos <- c(as.numeric(rownames(unique(data.frame(chrom = temp$chrom)[1]))), nrow(temp) + 1)
  numpos <- numeric(length(pos) - 1)
  for (i in 1:(length(pos) - 1)) {
    numpos[i] <- round((pos[i] + pos[i + 1] - 1) / 2)
  }
  
  mycols <- rep(NA, nrow(temp))
  mycols[1] <- colset[1]
  for (i in 2:nrow(temp)) {
    if (temp$chrom[i] == temp$chrom[i - 1]) {
      mycols[i] <- mycols[i - 1]
    } else {
      mycols[i] <- colset[which(colset != mycols[i - 1])]
    }
  }
  
  het_per_kb <- 1000 * temp[[hets_col]] / temp[[calls_col]]
  bp <- barplot(het_per_kb, col = mycols, border = mycols, main = plotname, ylim = c(0, 6))
  axis(1, at = bp[numpos], labels = FALSE)
  text(bp[numpos], par("usr")[3] - 0.5, labels = unique(temp$chrom), adj = 0.5, xpd = NA, cex = 0.8)
  mtext(side = 1, line = 3.5, "Chromosomes", font = 2)
  mtext(side = 2, line = 3.5, "Heterozygosity per kb", font = 2)
}

plotwinhethist <- function(samplename, labelx = FALSE, colset, temp) {
  calls_col <- paste0("calls_", samplename)
  hets_col  <- paste0("hets_",  samplename)
  mycol <- colset[1]
  
  h <- hist(1000 * temp[[hets_col]] / temp[[calls_col]], breaks = seq(0, 10, by = 0.1), plot = FALSE)
  bin_counts <- c(h$counts[1:59], sum(h$counts[60:length(h$counts)]))
  bp <- barplot(bin_counts, col = mycol, border = mycol, ylim = c(0, 300))
  
  if (labelx) {
    axis(1, at = seq(0, length(bp), by = 10), labels = 0:6)
  }
  mtext(side = 1, line = 3.5, "Het. per kb", font = 2)
  mtext(side = 2, line = 3.5, "# of windows", font = 2)
}

# ──────────────── FIRST PDF: All samples ────────────────
pdf("slidingwindowhets.pdf", width = 6.85, height = 5.5, pointsize = 8)
for (s in search_strings) {
  file_match <- all_files[grep(s, all_files)]
  if (!length(file_match)) { warning(paste("No file:", s)); next }
  data <- read.table(file_match[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("chrom","window_start","sites_total",
                      paste0("calls_", s), paste0("hets_", s))
  data$chrom <- trimws(data$chrom)
  data[[paste0("calls_", s)]] <- as.numeric(data[[paste0("calls_", s)]])
  data[[paste0("hets_", s)]]  <- as.numeric(data[[paste0("hets_", s)]])
  data <- na.omit(data)
  temp <- data[data[[paste0("calls_", s)]] >= 0.30 * winsize, ]
  if (!nrow(temp)) stop("No valid windows for ", s)
  
  if (grepl("ENP", s))       colset <- c("#1B9E77", "#A6D854")
  else if (grepl("GOC", s))  colset <- c("#D95F02", "#FC8D62")
  else if (grepl("ESP", s))  colset <- c("#7570B3", "#8DA0CB")
  else                       colset <- c("grey40","grey70")
  
  par(mfrow = c(2, 1), mar = c(4, 3, 2, 1))
  plotwinhet(s, colset, temp)
  plotwinhethist(s, labelx = TRUE, colset, temp)
}
dev.off()

# ──────────────── SECOND PDF: Three specific individuals ────────────────
selected_samples <- c("ESPCL10","ENPAK25","GOC002")
pdf("three_individuals.pdf", width = 10, height = 8, pointsize = 8)
par(mfrow = c(length(selected_samples), 2), mar = c(4, 3, 2, 1))
for (s in selected_samples) {
  file_match <- all_files[grep(s, all_files)]
  if (!length(file_match)) { warning(paste("No file:", s)); next }
  data <- read.table(file_match[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("chrom","window_start","sites_total",
                      paste0("calls_", s), paste0("hets_", s))
  data$chrom <- trimws(data$chrom)
  data[[paste0("calls_", s)]] <- as.numeric(data[[paste0("calls_", s)]])
  data[[paste0("hets_", s)]]  <- as.numeric(data[[paste0("hets_", s)]])
  data <- na.omit(data)
  temp <- data[data[[paste0("calls_", s)]] >= 0.30 * winsize, ]
  if (!nrow(temp)) { warning("No valid windows for ", s); next }
  
  if (grepl("ENP", s))       colset <- c("#1B9E77", "#A6D854")
  else if (grepl("GOC", s))  colset <- c("#D95F02", "#FC8D62")
  else if (grepl("ESP", s))  colset <- c("#7570B3", "#8DA0CB")
  else                       colset <- c("grey40","grey70")
  
  plotwinhet(s, colset, temp)
  plotwinhethist(s, labelx = TRUE, colset, temp)
}
dev.off()

# ──────────────── THIRD PDF: All sliding-window barplots only (3×4 layout, no clipping) ────────────────
skip_samples <- c("ENPOR12", "GOC077")
all_to_plot  <- setdiff(search_strings, skip_samples)
pdf("all_slidingwindows_3x4_fixed.pdf", width = 16, height = 10, pointsize = 8)
par(mfrow = c(3, 4), mar = c(5, 5, 3, 1), oma = c(0, 0, 0, 0))
for (s in all_to_plot) {
  file_match <- all_files[grep(s, all_files)]
  if (!length(file_match)) {
    warning("No file for sample:", s)
    next
  }
  data <- read.table(file_match[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("chrom", "window_start", "sites_total",
                      paste0("calls_", s), paste0("hets_", s))
  data$chrom <- trimws(data$chrom)
  data[[paste0("calls_", s)]] <- as.numeric(data[[paste0("calls_", s)]])
  data[[paste0("hets_", s)]]  <- as.numeric(data[[paste0("hets_", s)]])
  data <- na.omit(data)
  temp <- data[data[[paste0("calls_", s)]] >= 0.30 * winsize, ]
  if (!nrow(temp)) {
    warning("No valid windows for sample:", s)
    next
  }
  
  if (grepl("ENP", s))       colset <- c("#1B9E77", "#A6D854")
  else if (grepl("GOC", s))  colset <- c("#D95F02", "#FC8D62")
  else if (grepl("ESP", s))  colset <- c("#7570B3", "#8DA0CB")
  else                       colset <- c("grey40", "grey70")
  
  plotwinhet(s, colset, temp)
}
dev.off()

# ──────────────── FIFTH PDF: Three specific individuals (het-only barplots) ────────────────
selected_samples <- c("ESPCL10","ENPAK25","GOC002")
pdf("three_individuals_hetonly.pdf", width = 8, height = 12, pointsize = 10)
par(mfrow = c(length(selected_samples), 1), mar = c(5, 5, 3, 1))
for (s in selected_samples) {
  file_match <- all_files[grep(s, all_files)]
  if (!length(file_match)) {
    warning("No file for sample:", s)
    next
  }
  data <- read.table(file_match[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("chrom","window_start","sites_total",
                      paste0("calls_", s), paste0("hets_", s))
  data$chrom <- trimws(data$chrom)
  data[[paste0("calls_", s)]] <- as.numeric(data[[paste0("calls_", s)]])
  data[[paste0("hets_", s)]]  <- as.numeric(data[[paste0("hets_", s)]])
  data <- na.omit(data)
  temp <- data[data[[paste0("calls_", s)]] >= 0.30 * winsize, ]
  if (!nrow(temp)) {
    warning("No valid windows for sample:", s)
    next
  }
  
  if (grepl("ENP", s))       colset <- c("#1B9E77", "#A6D854")
  else if (grepl("GOC", s))  colset <- c("#D95F02", "#FC8D62")
  else if (grepl("ESP", s))  colset <- c("#7570B3", "#8DA0CB")
  else                       colset <- c("grey40","grey70")
  
  plotwinhet(s, colset, temp)
}
dev.off()

