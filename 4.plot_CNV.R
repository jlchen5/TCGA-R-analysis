##### Load Packages #####
library(tidyverse)
##### load data #####

dataCNV <- read.delim(file = "../TCGA-HNSC.gistic.tsv",header = T)

#### id convert  #####
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
dataCNV$Gene.Symbol <- as.character(dataCNV$Gene.Symbol)
dataCNV$Gene.Symbol <- sub("[.][0-9]*","",dataCNV$Gene.Symbol)

G_list2 <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=dataCNV$Gene.Symbol,mart= mart)
empty_rows <- G_list2[, 2] == ""
G_list2 <- G_list[!empty_rows, ]

dataCNV2 <- merge(dataCNV,G_list2,by.x="Gene.Symbol",by.y="ensembl_gene_id", no.dups = TRUE)
dataCNV2 <- dataCNV2[,c(528,1:527)]
rownames(dataCNV2) <- make.unique(as.character(dataCNV2[, 1]))

##### find PARG #####
f <- read.table("../PARG-overlap of PANoptosis and Tumor resistance genes.txt")	
dataCNV3 <- dataCNV2[dataCNV2$hgnc_symbol %in% f$V1, ]


########################
####### add anno #######
########################
group <- read.delim(file = '../TCGA-HNSC.GDC_sampletype.tsv',header = T)
rownames(group) <- group[,1]
rownames(group) <- gsub("-", ".", rownames(group))
#group <- group[,-1]
head(group)

dataCNV4 <- dataCNV3[,-c(1,2)]

########################
##### plot heatmap #####
########################
pheatmap::pheatmap(plotdata,
                   #scale='row',
                   annotation_col = group,
                   color = colorRampPalette(colors = c("blue","grey","red"))(3),
                   show_rownames=T,show_colnames = F,
                   cluster_cols = T,
                   cluster_rows = F)

