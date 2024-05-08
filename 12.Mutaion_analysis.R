# 载入R包
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(dplyr)
library(easyTCGA)
library(maftools)

# 读入数据
getsnvmaf("TCGA-HNSC")
load(file = "./output_snv/TCGA-HNSC_matched_maf_and_clin.rdata")



# cluster 1 and 2
load("ConsensusClusterPlus_results.rdata")
subtype_df <- data.frame(results[[2]]$consensusClass)
rownames(subtype_df) <- gsub(".",'-',rownames(subtype_df),fixed = T)
colnames(subtype_df) <- "subtype"
subtype_df$Sample_ID <- rownames(subtype_df)

subtype_df$Tumor_Sample_Barcode <- gsub("-01[AB]", "", subtype_df$Sample_ID)
subtype_df$Tumor_Sample_Barcode2 <- gsub("-11A", "", subtype_df$Tumor_Sample_Barcode)
subtype_df$Tumor_Sample_Barcode3 <- gsub("-06A", "", subtype_df$Tumor_Sample_Barcode2)

subtype_df2 <- subtype_df[,c(1,5)]
colnames(subtype_df2) <- c("subtype","barcode")
  
clin_snv_anno <- merge(clin_snv,subtype_df2,
                       by.x = "submitter_id",
                       by.y="barcode")
clin_snv_anno <- clin_snv_anno[,c(2,1,3:76)]



# 直接读取
coad_maf <- read.maf(data,clin_snv_anno,isTCGA = T)

maftools::oncoplot(coad_maf,clinicalFeatures = "subtype",
                   sortByAnnotation = TRUE)

somaticInteractions(coad_maf, pvalue = c(0.05, 0.01, 0.001),
                    top=20,fontSize = 0.5)
