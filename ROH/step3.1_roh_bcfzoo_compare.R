
# Adjusted the previous script to actually plot the line equations neatly. Same as original script only plotting changes. 

rm(list=ls()); cat("\014"); options(echo=TRUE)
setwd("/data/shared/snigenda/finwhale_projects/fin_genomics/ROH/")
library(dplyr); library(ggplot2); library(reshape2); library(tidyr)

get_rohnames <- function(rohlens) {
  out=character(0)
  for(i in seq_along(rohlens)) {
    tmp <- if(i==length(rohlens)) paste0(rohlens[i],"_Inf") else paste0(rohlens[i],"_",rohlens[i+1])
    out <- c(out,tmp)
  }
  out
}

get_forplot <- function(rohsum) {
  fp <- rohsum %>%
    select(-starts_with("froh"),-GenomeHet) %>%
    melt(id.vars=c("SampleId","PopId"))
  dt <- colsplit(fp$variable,pattern="_",names=c("type","software","length"))
  fp <- cbind(fp,dt)
  b <- fp %>% filter(software=="bcf") %>% rename(bcftools=value) %>% select(SampleId,PopId,type,length,bcftools)
  z <- fp %>% filter(software=="zoo") %>% rename(RZooRoH=value) %>% select(SampleId,PopId,type,length,RZooRoH)
  full_join(b,z,by=c("SampleId","PopId","type","length")) %>%
    drop_na() %>%
    mutate(type=ifelse(type=="N","Total number","Total length (Mb)"))
}

plotdir <- "./plots/"; outdir <- "./derived_data/"
dir.create(plotdir,recursive=TRUE,showWarnings=FALSE)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
today <- format(Sys.Date(),"%Y%m%d")
genomelen <- 2239549461
rohlens <- c(0.1,1,5)*1e+6
rohcats <- get_rohnames(c(0.1,1,5))
sessionInfo()
source("/data/shared/snigenda/finwhale_projects/scripts/winfingenomics/config/plotting_config.R")
rohsum2 <- readRDS("./derived_data/rohsummary_zoobcf_20250603.rds")

forplot <- get_forplot(rohsum2)
types <- unique(forplot$type)
lens  <- unique(forplot$length)
lenslab <- c("[0.1, 1) Mb","[1, 5) Mb","[5, Inf) Mb","[0.1, Inf) Mb")
lenscat <- c("short","medium","long","all")

plotlist <- vector("list",length(types)*length(lens)); counter <- 0

for(ii in types){
  for(j in seq_along(lens)){
    jj <- lens[j]; lenlab <- lenslab[j]; lencat <- lenscat[j]; counter <- counter+1
    dat <- forplot %>% filter(type==ii, length==jj)
    if(ii=="Total length (Mb)") dat <- dat %>% mutate(bcftools=bcftools/1e+6, RZooRoH=RZooRoH/1e+6)
    pr <- range(dat$bcftools,dat$RZooRoH); dr <- diff(pr); mgn <- 0.05*dr

    mall <- lm(RZooRoH~bcftools,data=dat)
    a1 <- coef(mall)[1]; b1 <- coef(mall)[2]; r1 <- summary(mall)$r.squared
    lab_all <- paste0("italic(y)==",round(a1,2)," + ",round(b1,2)," %.% italic(x)~~italic(R)^2==",round(r1,2))
    xo <- pr[2]-mgn; yo <- pr[1]+mgn
    df_all <- data.frame(x=xo,y=yo,label=lab_all)

    pops <- sort(unique(dat$PopId)); np <- length(pops)
    xp <- pr[1] + 0.005*dr
    off <- 0.04*dr
    raise <- 0.02*dr
    yv <- (pr[2]-mgn + raise) - (seq_along(pops)-1)*off
    df_pe <- data.frame(PopId=pops,x=xp,y=yv,label=NA_character_)
    for(k in seq_along(pops)){
      pp <- pops[k]
      sub <- dat %>% filter(PopId==pp)
      mpp <- lm(RZooRoH~bcftools,data=sub)
      a2 <- coef(mpp)[1]; b2 <- coef(mpp)[2]; r2 <- summary(mpp)$r.squared
      df_pe$label[k] <- paste0("italic(y)==",round(a2,2)," + ",round(b2,2)," %.% italic(x)~~italic(R)^2==",round(r2,2))
    }

    fml <- y~x
    pp <- ggplot() +
      geom_point(data=dat,aes(x=bcftools,y=RZooRoH),shape=".",color="gray60") +
      geom_smooth(data=dat,aes(x=bcftools,y=RZooRoH),method="lm",se=FALSE,formula=fml,size=0.5,linetype="dashed",color="black") +
      geom_smooth(data=dat,aes(x=bcftools,y=RZooRoH,color=PopId),method="lm",se=FALSE,formula=fml,size=0.8) +
      geom_point(data=dat,aes(x=bcftools,y=RZooRoH,color=PopId)) +
      geom_text(data=df_all,aes(x=x,y=y,label=label),parse=TRUE,hjust=1,vjust=0,color="black",size=3) +
      geom_text(data=df_pe,aes(x=x,y=y,label=label,color=PopId),parse=TRUE,hjust=0,vjust=1,size=3) +
      geom_abline(slope=1,intercept=0,color="darkgray",linetype="dotted") +
      labs(title=paste(ii,"of",lencat,"ROH"),subtitle=paste("ROH length:",lenlab)) +
      coord_fixed(ratio=1,xlim=pr,ylim=pr) +
      scale_color_manual(values=mycolors) +
      theme_bw() +
      theme(legend.position="none",plot.title=element_text(size=12),plot.subtitle=element_text(size=10))

    plotlist[[counter]] <- pp
  }
}

ppout <- ggpubr::ggarrange(plotlist=plotlist,nrow=length(types),ncol=length(lens))
ggsave(paste0("FigureS8.ROH_bcfzoo_compare_",today,".pdf"),path=plotdir,width=14,height=8)
write.csv(forplot,file="./derived_data/FigS9.csv")
date(); closeAllConnections()
