##### rawdata Analysis #####
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(ggpubr)

##### load data #####
rawdata <- read.delim(file = "../TCGA-HNSC.htseq_counts.tsv",header = T)

#ID convert from ensembl_gene_id to gene_symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
rawdata$Ensembl_ID <- as.character(rawdata$Ensembl_ID)
rawdata$Ensembl_ID <- sub("[.][0-9]*","",rawdata$Ensembl_ID)

G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=rawdata$Ensembl_ID,mart= mart)
empty_rows <- G_list[, 2] == ""
G_list <- G_list[!empty_rows, ]

rawdata2 <- merge(rawdata,G_list,by.x="Ensembl_ID",by.y="ensembl_gene_id", no.dups = TRUE)

count <- rawdata2[,c(548,1:547)]
rownames(count) <- make.unique(as.character(count[, 1]))
count <- count[,-c(1,2)]


##### get pheno #####
pheno.raw <- read.delim("../TCGA-HNSC.GDC_sampletype.tsv",header = T)
rownames(pheno.raw) <- pheno.raw[,1]
rownames(pheno.raw) <- gsub("-", ".", rownames(pheno.raw))
pheno.raw$submitter_id.samples <- gsub("-", ".", pheno.raw$submitter_id.samples)
pheno <- pheno.raw[pheno.raw$submitter_id.samples %in% colnames(count),] 
# rownames(pheno) <- pheno$pheno.raw$submitter_id.samples
# count <- count[,pheno$pheno.raw$submitter_id.samples]
# pheno  <- pheno[colnames(count),]
pheno

table(colnames(count)==pheno$submitter_id.samples)
colnames(count) <- pheno$submitter_id.samples

##### PCA #####
colData <- data.frame(row.names =colnames(count),type=pheno$sample_type.samples)
dds <- DESeqDataSetFromMatrix(countData = round(count), colData = colData, design = ~ type) 
cpm <- fpm(dds,robust = FALSE)
cpm <- as.data.frame(cpm)
cpm$gene_name <- rownames(cpm)

save(cpm,file = "cpm.rdata")

vst <- vst(dds)
plotPCA(vst,intgroup="type") 
data <- plotPCA(vst,intgroup="type",returnData = TRUE)
data$ID <- rownames(data)
data.pheno <- merge(data,pheno,by=0)
data.pheno

ggplot(data.pheno, aes(x = PC1, y=PC2)) +
  geom_point(size=2,stroke = 1,aes(color=group)) + 
  labs(title=" ") + theme_bw()+
  theme(plot.title=element_text(hjust=0.5),title =element_text(size=12) )  +
  xlab(paste0( " PC1: " , " 19% variance " ))+
  ylab(paste0( " PC2: " , " 14% variance " ))



dds <- DESeq(dds)
res <- results(dds,contrast = c("type","Primary Tumor","Solid Tissue Normal"))
res <- as.data.frame(res)
res$gene_name <- rownames(res)
save(res,file = 'deg_res.rdata')

cpm.res <- as.data.frame(merge(cpm,res,by="gene_name"))
cpm.res$gene_name <- rownames(res)
#cpm.res.anno <- merge(anno,cpm.res,by="repeatMasker.repName")
write.csv(cpm.res,file="HNSC_DEG.csv")

f <- read.table("../PARG-overlap of PANoptosis and Tumor resistance genes.txt")	
res_20 <- res[res$gene_name %in% f$V1, ]


# 创建正负柱状图
barplot(res_20$log2FoldChange, names.arg = res_20$gene_name ,las=2,
        col = ifelse(res_20$log2FoldChange >= 0, "blue", "red"),
        main = "PARG Genes", ylab = "Log2(FoldChange)")


##### correlation #####
# calculate correlation

data_corr <- cor(t(plotdata))
data_corr <- as.data.frame(data_corr)


pheatmap::pheatmap(data_corr,cluster_cols  = T,fontsize_row=12,
         color = colorRampPalette(c("#4575b4","#abd9e9","#abd9e9","#ffffbf","#fee090","#fee090","#d73027"))(100),
         fontsize_col = 12,angle_col = 315,main = "PARG Correlation")






