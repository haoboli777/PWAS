rm(list=ls())
setwd("E:/projects/Project/")
data = data.table::fread("./2.gwas-data/ukb_new_gwas.txt",header = T,fill=T)
head(data)
data = data[,c("RSID","CHR","BP","A1","A0","BETA","SE","P","NMISS")]
colnames(data) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N")
write.table(data,"./2.gwas-data/gwas-clean",sep = "\t",quote = F,row.names = F)
