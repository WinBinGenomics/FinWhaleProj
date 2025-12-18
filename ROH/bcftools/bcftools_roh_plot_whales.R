## This script is for filtering and plotting bcftools roh output
## Input is (optionally gzipped) .out files from bcftools ROH
## subsetted down to have one file for each individual as the total file size can be quite large

# Author: Chris Kyriazis, modified by Sergio Nigenda for fin whale data
# Adjusted by Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale project.

#### Set working directory, file names, individuals, etc ####

#setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/GOC")
#setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/ENP")
setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/ESP")

# set filename excluding individual and ".out.gz"
# file <- "_GOC_concat_fwhale_roh_bcftools_G30_ACANGT"
#file <- "_ENP_concat_fwhale_roh_bcftools_G30_ACANGT"
file <- "_ESP_concat_fwhale_roh_bcftools_G30_ACANGT"

## set array of individuals 
# individuals <- c("GOC002","GOC006","GOC025","GOC038","GOC050","GOC053","GOC063","GOC068","GOC071","GOC077","GOC080","GOC082","GOC086","GOC091","GOC100","GOC111","GOC112","GOC116","GOC125")
#individuals <- c("ENPAK19","ENPAK20","ENPAK21","ENPAK22","ENPAK23","ENPAK24","ENPAK25","ENPAK26","ENPAK27","ENPAK28","ENPAK29","ENPAK30","ENPBC16","ENPBC17","ENPBC18","ENPCA01","ENPCA02","ENPCA03","ENPCA04","ENPCA05","ENPCA06","ENPCA07","ENPCA08","ENPOR10","ENPOR11","ENPOR13","ENPWA14","ENPWA15")
individuals <- c("ESPCL01","ESPCL02","ESPCL03","ESPCL04","ESPCL05","ESPCL06","ESPCL07","ESPCL08","ESPCL09","ESPCL10","ESPCL11","ESPCL12","ESPCL13","ESPCL14","ESPCL15","ESPCL16","ESPCL17","ESPCL18","ESPCL19","ESPCL20")

## need genome length as denominator for Froh
## might alternatively consider using total # callable sites

# Minke total length of scaffolds larger than 1Mb used for analyses
#genome_length <- 2324.429847

#blue whale autosome length
genome_length <- 2239.549461 


# min size allowable for ROHs - typically 100kb or 300kb
min_roh_length=100000


## Define function that takes in data frame and divides it into three length classes
classify_roh <- function(roh_dataframe, min_roh_length){
  short_roh <- subset(roh_dataframe,length>min_roh_length & length<1000000) # roh_dataframe[length>100000 & length<1000000]
  med_roh <- subset(roh_dataframe, length>1000000 & length<10000000)
  long_roh <-  subset(roh_dataframe, length>10000000 & length<100000000)
  
  #sum each class and divide by 1000000 to convert to Mb
  sum_short_Mb <- sum(short_roh$length)/1000000
  sum_med_Mb <- sum(med_roh$length)/1000000
  sum_long_Mb <- sum(long_roh$length)/1000000
  
  print(paste("This individual has",dim(short_roh)[1],"short ROHs summing to",sum_short_Mb, "Mb",
              dim(med_roh)[1],"medium ROHs summing to",sum_med_Mb,"Mb, and",
              dim(long_roh)[1],"long ROHs summing to",sum_long_Mb, "Mb"))
  
  return(c(sum_short_Mb, sum_med_Mb, sum_long_Mb))
  #roh_matrix_Mb <- roh_matrix/1000000 #convert to Mb
  
}

## Define function to read in output files for each individual and filter out ROHs less than min_roh_length
read_filter_roh <- function(data, min_roh_length){
  output <- read.table(paste(data,".out.gz",sep=""), col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"), fill=T)
  output1 <- subset(output, row_type == "RG")
  output_class_sums <- classify_roh(output1, min_roh_length=min_roh_length)
  return(output_class_sums)
}



#### Read in data, filter, and classify ROHs ####

## initialize data frame
roh_size_df <- data.frame(matrix(nrow=3, ncol=length(individuals)))
colnames(roh_size_df) <- individuals
froh <- c()

## read in data for each individaul
## note that this can be VERY slow - takes several min per individual
for(i in 1:length(individuals)){
  roh_size_df[,i] <- read_filter_roh(paste(individuals[i],file, sep=""), min_roh_length=min_roh_length)
  froh = c(froh,sum(roh_size_df[2:3,i])/genome_length) # sum ROHs > 1Mb and divide by genome length to estimate Froh
}

froh

#### Plot results ####

setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/rohbcftools/plots")

png(paste("barplot",file,"_subset.png",sep=""), width=16, height=10, units="in", res=600)

barplot(as.matrix(roh_size_df), names.arg = individuals, col=c("seashell","seashell3", "seashell4"), ylab="Summed ROH length (Mb)", xlab="Individual", ylim=c(0,1000))
legend("topright",legend=c(paste(min_roh_length/1000000,"-1 Mb",sep=""), "1-10 Mb", "10-100 Mb"), col=c("red","orange", "yellow"), fill=c("seashell","seashell3", "seashell4"))

dev.off()
