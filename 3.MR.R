### MR分析
rm(list=ls())
library(dplyr)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)
setwd("E:/projects/Project/")

######################
# 获取工具变量
######################
rm(list=ls())
gwas = data.table::fread("./2.gwas-data/gwas-clean",header = T)
gwas = gwas[,c("SNP","A1","A2","BETA","SE","P")]
pwas = data.table::fread("./4.res-pwas/res-pwas",header = T)
pwas = pwas[pwas$FDR==TRUE,]
pwas = unique(pwas$FILE)
for (l in 1:length(pwas)){
  #读取暴露
  pqtl = data.table::fread(sprintf("./3.pqtl-data/EA/%s.PHENO1.glm.linear",pwas[l]),header = T)
  pqtl1 = pqtl[pqtl$REF == pqtl$A1,]
  pqtl2 = pqtl[pqtl$ALT == pqtl$A1,]
  pqtl1 = pqtl1[,c("ID","#CHROM","POS","REF","ALT","BETA","SE","P","OBS_CT","A1_FREQ")]
  colnames(pqtl1) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","MAF")
  pqtl2 = pqtl2[,c("ID","#CHROM","POS","ALT","REF","BETA","SE","P","OBS_CT","A1_FREQ")]
  colnames(pqtl2) = c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","MAF")
  pqtl = rbind(pqtl1,pqtl2)
  # 数据格式化
  exp_dat = format_data(pqtl, snp_col = "SNP",
                        beta_col = "BETA",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        eaf_col = "MAF")
  exp_dat = try(clump_data(exp_dat,  clump_kb =  1000,clump_r2 =  0.1,clump_p1 =  5e-5,pop =  "EUR" ))
  exp_dat = as.data.frame(exp_dat[,1])
  colnames(exp_dat) = "SNP"
  exp_dat = merge(exp_dat,pqtl,by="SNP")
  # 对齐结局
  x_y = merge(exp_dat,gwas,by="SNP")
  x_y = x_y %>% filter((A1.x==A1.y & A2.x==A2.y)|(A1.x==A2.y & A2.x==A1.y)) 
  x_y = x_y %>% mutate(BETA.y = ifelse(A1.x==A1.y,BETA.y,-BETA.y))
  
  x_y = x_y[,c("SNP","CHR","BP","N","MAF","A1.x","A2.x","BETA.x","SE.x","P.x","BETA.y","SE.y","P.y")]

  write.table(x_y,sprintf("./5.res-mr/data-iv-5e5-0.1/%s",pwas[l]),row.names = F,quote = F,sep = "\t")
}


######################
# 进行MR分析
######################
rm(list=ls())
# 读取数据
pwas = data.table::fread("./4.res-pwas/res-pwas",header = T)
pwas = pwas[pwas$FDR==TRUE,]

# 创建结果data.frame
result1 <- data.frame(matrix(NA,length(pwas$FILE),14))
colnames(result1) <- c("method","beta","se","CILower","CIUpper","OR","ORLower","ORUpper","pval","seq","protein","num_SNP","F","PVE")
result2 <- data.frame(matrix(NA,length(pwas$FILE),14))
colnames(result2) <- c("method","beta","se","CILower","CIUpper","OR","ORLower","ORUpper","pval","seq","protein","num_SNP","F","PVE")
result3 <- data.frame(matrix(NA,length(pwas$FILE),14))
colnames(result3) <- c("method","beta","se","CILower","CIUpper","OR","ORLower","ORUpper","pval","seq","protein","num_SNP","F","PVE")
n=4

# for循环执行蛋白质和疾病的MR分析
for (l in 1:length(pwas$FILE)){
  print(l)
  
  # 读取蛋白质工具变量
  data = data.table::fread(sprintf("./5.res-mr/data-iv-5e5-0.1/%s",pwas[l,1]),header = T)

  # 计算R2和F statstic
  data$pve = 2*(data$BETA.x^2)*data$MAF*(1-data$MAF)
  R2 = sum(data$pve)
  N <- max(data$N) 
  k <- length(data$SNP)
  Fstat <- R2*(N-k-1)/((1-R2)*k) 
  print(c(pwas[l,1],R2,Fstat,length(data$SNP)))
  
  mr_input <- mr_input(bx = data$BETA.x, bxse = data$SE.x,by = data$BETA.y, byse = data$SE.y,snps = data$SNP)
  
  # MR分析 IVW算法
  if (length(data$SNP)>0){
    mr_ivw = mr_ivw(mr_input)
    result1[l,1] <- "Inverse-variance weighted"
    result1[l,2] <- round(mr_ivw@Estimate,n)
    result1[l,3] <- round(mr_ivw@StdError,n)
    result1[l,4] <- round(mr_ivw@CILower,n)
    result1[l,5] <- round(mr_ivw@CIUpper,n)
    result1[l,6] <- round(exp(mr_ivw@Estimate),n)
    result1[l,7] <- round(exp(mr_ivw@CILower),n)
    result1[l,8] <- round(exp(mr_ivw@CIUpper),n)
    result1[l,9] <- mr_ivw@Pvalue
    result1[l,10] <- pwas[l,1]
    result1[l,11] <- pwas[l,2]
    result1[l,12] <- length(data$SNP)
    result1[l,13] <- Fstat
    result1[l,14] <- R2
  }else{
    result1[l,10] <- pwas[l,1]
    result1[l,11] <- pwas[l,2]
  }
  
  # MR分析 MR egger算法
  if (length(data$SNP)>2){
    MR_Egger = mr_egger(mr_input)
    result2[l,1] <- "EGGER"
    result2[l,2] <- round(MR_Egger@Estimate,n)
    result2[l,3] <- round(MR_Egger@StdError.Est,n)
    result2[l,4] <- round(MR_Egger@CILower.Est,n)
    result2[l,5] <- round(MR_Egger@CIUpper.Est,n)
    result2[l,6] <- round(exp(MR_Egger@Estimate),n)
    result2[l,7] <- round(exp(MR_Egger@CILower.Est),n)
    result2[l,8] <- round(exp(MR_Egger@CIUpper.Est),n)
    result2[l,9] <- MR_Egger@Pvalue.Est
    result2[l,10] <- pwas[l,1]
    result2[l,11] <- pwas[l,2]
    result2[l,12] <- length(data$SNP)
    result2[l,13] <- Fstat
    result2[l,14] <- R2
    
    result3[l,1] <- "INT"
    result3[l,2] <- round(MR_Egger@Intercept,n)
    result3[l,3] <- round(MR_Egger@StdError.Int,n)
    result3[l,4] <- round(MR_Egger@CILower.Int,n)
    result3[l,5] <- round(MR_Egger@CIUpper.Int,n)
    result3[l,6] <- round(exp(MR_Egger@Intercept),n)
    result3[l,7] <- round(exp(MR_Egger@CILower.Int),n)
    result3[l,8] <- round(exp(MR_Egger@CIUpper.Int),n)
    result3[l,9] <- MR_Egger@Pvalue.Int
    result3[l,10] <- pwas[l,1]
    result3[l,11] <- pwas[l,2]
    result3[l,12] <- length(data$SNP)
    result3[l,13] <- Fstat
    result3[l,14] <- R2
  }else{
    result2[l,1] <- "EGGER"
    result2[l,11] <- pwas[l,1]
    result3[l,1] <- "INT"
    result3[l,11] <- pwas[l,1]
  }
}
write.table(result1,"./5.res-mr/5E5-0.1-IVW",sep = "\t",row.names = F,quote = F)
write.table(result2,"./5.res-mr/5E5-0.1-MRE",sep = "\t",row.names = F,quote = F)
write.table(result3,"./5.res-mr/5E5-0.1-MRE-INT",sep = "\t",row.names = F,quote = F)
