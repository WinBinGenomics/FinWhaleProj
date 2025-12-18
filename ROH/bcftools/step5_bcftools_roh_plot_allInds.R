# Bar plots to compare ROH between the populations + calculate FOH values. 

# Authors: Kaden Winspear and Sergio Nigenda @ CSUSM. 


# parent directory containing the three pop folders
setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools")

file <- "_concat_fwhale_roh_bcftools_G30_ACANGT"

## ALL individuals, named "ID_POP"
individuals <- c(
  # ENP #########################
  "ENPAK19_ENP","ENPAK20_ENP","ENPAK21_ENP","ENPAK22_ENP","ENPAK23_ENP","ENPAK24_ENP",
  "ENPAK25_ENP","ENPAK26_ENP","ENPAK27_ENP","ENPAK28_ENP","ENPAK29_ENP","ENPAK30_ENP",
  "ENPBC16_ENP","ENPBC17_ENP","ENPBC18_ENP","ENPCA02_ENP","ENPCA03_ENP","ENPCA04_ENP",
  "ENPCA05_ENP","ENPCA06_ENP","ENPCA07_ENP","ENPCA08_ENP","ENPOR10_ENP","ENPOR11_ENP",
  "ENPOR13_ENP","ENPWA14_ENP","ENPWA15_ENP",
  # ESP ##########################
  "ESPCL01_ESP","ESPCL02_ESP","ESPCL03_ESP","ESPCL04_ESP","ESPCL05_ESP","ESPCL06_ESP",
  "ESPCL07_ESP","ESPCL08_ESP","ESPCL09_ESP","ESPCL10_ESP","ESPCL11_ESP","ESPCL12_ESP",
  "ESPCL13_ESP","ESPCL14_ESP","ESPCL15_ESP","ESPCL16_ESP","ESPCL17_ESP","ESPCL18_ESP",
  "ESPCL19_ESP","ESPCL20_ESP",
  # GOC ##########################
  "GOC002_GOC","GOC006_GOC","GOC025_GOC","GOC038_GOC","GOC050_GOC","GOC053_GOC",
  "GOC063_GOC","GOC068_GOC","GOC071_GOC","GOC080_GOC","GOC082_GOC","GOC086_GOC",
  "GOC091_GOC","GOC100_GOC","GOC111_GOC","GOC112_GOC","GOC116_GOC","GOC125_GOC"
)
# autosomes length
genome_length   <- 2239.549461
min_roh_length <- 100000

#### Functions (with progress reporting) ####
classify_roh <- function(roh_dataframe, min_roh_length){
  short_roh <- subset(roh_dataframe, length > min_roh_length & length <= 1e6)
  med_roh   <- subset(roh_dataframe, length > 1e6           & length <= 5e6)
  long_roh  <- subset(roh_dataframe, length > 5e6)

  sum_short_Mb <- sum(short_roh$length)/1e6
  sum_med_Mb   <- sum(med_roh$length)  /1e6
  sum_long_Mb  <- sum(long_roh$length) /1e6

  # progress report
  print(paste("This individual has", dim(short_roh)[1], "short ROHs summing to", sum_short_Mb, "Mb;",
              dim(med_roh)[1],   "medium ROHs summing to",   sum_med_Mb,   "Mb;",
              dim(long_roh)[1],  "long ROHs summing to",     sum_long_Mb,  "Mb"))

  return(c(sum_short_Mb, sum_med_Mb, sum_long_Mb))
}

read_filter_roh <- function(data, min_roh_length){
  output <- read.table(paste0(data, ".out.gz"),
                       col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"),
                       fill=TRUE)
  output1 <- subset(output, row_type == "RG")
  classify_roh(output1, min_roh_length)
}

#### Read, classify, calculate FROH ####

roh_size_df <- data.frame(matrix(nrow=3, ncol=length(individuals)))
colnames(roh_size_df) <- individuals
froh <- numeric(length(individuals))

for(i in seq_along(individuals)){
  parts <- strsplit(individuals[i], "_")[[1]]
  base  <- parts[1]
  pop   <- parts[2]
  path  <- file.path(pop, paste0(base, "_", pop, file))

  roh_size_df[,i] <- read_filter_roh(path, min_roh_length)
  froh[i]         <- sum(roh_size_df[2:3,i]) / genome_length
}
names(froh) <- individuals

# Order Populations
pop_order <- c("ESP","ENP","GOC")
order_idx <- unlist(lapply(pop_order, function(p) grep(paste0("_",p,"$"), colnames(roh_size_df))))
roh_size_df <- roh_size_df[, order_idx]
froh       <- froh[order_idx]
names(froh) <- colnames(roh_size_df)

# write combined FROH
write.table(data.frame(Individual=sub("_(ENP|ESP|GOC)$", "", names(froh)),
                       FROH=round(froh,5)),
            file="FROH_all_pops.txt",
            sep="\t", row.names=FALSE, quote=FALSE)

#### Plotting 
dir.create("plots", showWarnings=FALSE)
png("plots/ROH_barplot_all_pops.png", width=16, height=12, units="in", res=600)

# adjust margins
par(mar=c(12, 6, 4, 2))

xlabs <- sub("_(ENP|ESP|GOC)$", "", colnames(roh_size_df))

# Colors
zissou_cols <- hcl.colors(5, "Zissou 1")
cols <- c(zissou_cols[3],  # yellow (small runs)
          zissou_cols[4],  # orange (medium runs)
          zissou_cols[5])  # red (long runs)

# draw barplot and save midpoints for x-axis label
bp <- barplot(as.matrix(roh_size_df),
              horiz	 = FALSE,
              names.arg  = xlabs,
              las        = 2,        # rotate x labels vertical
              col        = cols,
              border     = NA,
              ylab	 = "Summed ROH length (Mb)",
              cex.names  = 0.6,
              cex.axis   = 1.2,
              cex.lab    = 1.4,
              font.axis  = 2,
              font.lab   = 2,
              ylim	 = c(0, max(colSums(roh_size_df))*1.05)
)

# add custom X-axis label below tick labels
mtext("Individuals", side = 1, line = 10, cex = 1.5, font = 2)

legend(x = "topleft", inset = c(0.02, 0), xpd = TRUE,
       legend    = c("0.1–1 Mb","1–5 Mb",">5 Mb"),
       fill      = cols,
       border    = NA,
       cex       = 1.5,
       text.font = 2,
       box.lwd   = 2
)

dev.off()
