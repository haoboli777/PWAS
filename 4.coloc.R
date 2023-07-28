#共定位分析
rm(list=ls())
library(dplyr)
library(coloc)
setwd("E:/projects/Project/")

# 读取gwas数据和PWAS结果
gwas = data.table::fread("./2.gwas-data/gwas-clean",header = T)
gwas = gwas[,c("SNP","A1","A2","BETA","SE","P")]
pwas = data.table::fread("./4.res-pwas/res-pwas",header = T)
pwas = pwas[pwas$FDR==TRUE,]

# GWAS的样本量和case比例
#N = 10000
#S = 0.3

# 创建结果data.frame
colocsum = data.frame()
colocres = data.frame()

# for循环执行蛋白质和疾病的coloc分析
for (l in 1:length(pwas$FILE)){
  print(l)
  
  # 读取pqtl数据
  pqtl = data.table::fread(sprintf("./3.pqtl-data/EA/%s.PHENO1.glm.linear",as.character(pwas[l,1])),header = T)
  pqtl1 = pqtl[pqtl$REF == pqtl$A1,]
  pqtl2 = pqtl[pqtl$ALT == pqtl$A1,]
  pqtl1 = pqtl1[,c("ID","#CHROM","POS","REF","ALT","BETA","SE","P","OBS_CT","A1_FREQ")]
  colnames(pqtl1) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","MAF")
  pqtl2 = pqtl2[,c("ID","#CHROM","POS","ALT","REF","BETA","SE","P","OBS_CT","A1_FREQ")]
  colnames(pqtl2) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","MAF")
  pqtl = rbind(pqtl1,pqtl2)
  
  # 对齐效应等位基因
  data = merge(pqtl,gwas,by = "SNP")
  data = data %>% filter((A1.x==A1.y & A2.x==A2.y)|(A1.x==A2.y & A2.x==A1.y)) 
  data = data %>% mutate(BETA.y = ifelse(A1.x==A1.y,BETA.y,-BETA.y))
  data = data[!duplicated(data$SNP),]
  
  data$var.x = data$SE.x**2
  data$var.y = data$SE.y**2
  data = data[data$var.x!=0,]
  data = data[data$var.y!=0,]
  
  write.table(data,sprintf("./6.res-coloc/data/%s",pwas[l,1]),row.names = F,quote = F,sep = "\t")
  
  # coloc要求数据准备

  data1 = data[,c("BETA.x","var.x","SNP","MAF","N")]
  data2 = data[,c("BETA.y","var.y","SNP")]
  
  colnames(data1) = c("beta","varbeta","snp","MAF","N")
  colnames(data2) = c("beta","varbeta","snp")
  
  data1 = as.list(data1)
  data2 = as.list(data2)
  #data2$N = N
  #data2$s = S
  
  data1$type = "quant"
  data2$type = "cc"
  
  # 共定位分析
  coloc = coloc.abf(data1,data2,p1=1e-4,p2=1e-4,p12=1e-5)

  
  # 结果汇总
  coloc_sum = as.data.frame(t(coloc$summary))
  coloc_res = coloc$results
  
  coloc_sum$gene = as.character(pwas[l,1])
  coloc_res$gene = as.character(pwas[l,1])
  
  colocsum = rbind(colocsum,coloc_sum)
  colocres = rbind(colocres,coloc_res)
  
}
write.table(colocsum,"./6.res-coloc/COLOC-SUM",sep = "\t",row.names = F,quote = F)
write.table(colocres,"./6.res-coloc/COLOC-RES",sep = "\t",row.names = F,quote = F)
