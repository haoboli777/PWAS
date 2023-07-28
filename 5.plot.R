
################################
### 绘图
################################
# manhattan plot
rm(list=ls())
library("CMplot")
library("stringr")
setwd("D:/project5/7.plot/manhattan/")
data = read.table("./PWAS1",sep = "\t",header=TRUE)
trait = unique(data$Trait)
i=1
for (i in 1:length(trait)){
  dat = data[data$Trait==trait[i],]
  dat1  = dat[,c("GENE","CHR","P0","PWAS.P")]
  dat2  = dat[dat$SIG==TRUE,]$GENE
  CMplot(dat1, plot.type="m", LOG10=TRUE, ylim=NULL, band = 0.2,
         # plot size
         bin.size=1e6,cex = 0.5,amplify=T,signal.cex =0.5,signal.col="red",
         # X/Y lable size
         cex.axis = 1.5,cex.lab=1.5,ylab="-log10(PWAS.P)",
         # highlight
         highlight=dat2, highlight.type = "p", highlight.cex=0.8, highlight.col="black",
         #highlight.text=dat2, highlight.text.cex=1.2, highlight.text.col="black",
         chr.den.col=NULL,
         file="pdf",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6,
         memo=sprintf("no-%s",trait[i]),main=sprintf("no-%s",trait[i]))
}

# ggforestplot
#devtools::install_github("NightingaleHealth/ggforestplot")
library(tidyverse)
library(ggforestplot)
rm(list=ls())
setwd("D:/project5/7.plot/forestplot/")

data = read.table("./data-1.txt",header = T,sep = "\t")
trait = unique(data$Trait)
for (i in 1:length(trait)){
  dat = data[data$Trait==trait[i],]
  p = forestplot(df = dat, name = TG, colour = class,estimate = beta, logodds = TRUE, title = sprintf("MR-%s",trait[i]), xlab = "OR")
  p = p+ theme(text=element_text(size=16,  family="serif"))
  pdf(sprintf("./%s.pdf",trait[i]),width = 8, height = 1+length(dat$TP)*0.2)
  print(p)
  dev.off()
}

# locuscomparer
rm(list=ls())
library(locuscomparer)
library(reshape)
?locuscompare
setwd("D:/project5/7.plot/locuscom/")
data = read.table("./data.txt",header = T,sep = "\t")
i=1
for (i in 1:length(data$TP)){
  
  locus = read.table(sprintf("../../6.res-coloc/data-P2/%s..%s",data[i,2],data[i,3]), header=T)
  
  gwas1 = locus[,c("SNP","P.x")]
  colnames(gwas1) = c("rsid","pval")
  
  gwas2 = locus[,c("SNP","P.y")]
  colnames(gwas2) = c("rsid","pval")

  
  snp = as.character(data[i,10])
  
  RES = locuscompare(in_fn1=gwas1, in_fn2=gwas2, snp=snp, population="EUR",combine=F,
                     title1="-log(P.pqtl)",title2="-log(P.gwas)",legend_position="topleft")
  RES$locuscompare
  
  pdf(sprintf("./B-%s.%s.pdf",data[i,2],data[i,3]))
  print(RES$locuscompare)
  dev.off()

}