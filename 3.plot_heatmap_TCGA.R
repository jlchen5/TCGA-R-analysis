########################
#### load packages #####
########################
rm(list = ls())
library('biomaRt')
library(dplyr)

########################
###### read data #######
########################
deg <- read.table(file = '../TCGA-HNSC.htseq_fpkm.tsv',header = T)

#ID convert from ensembl_gene_id to gene_symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
deg$Ensembl_ID <- as.character(deg$Ensembl_ID)
deg$Ensembl_ID <- sub("[.][0-9]*","",deg$Ensembl_ID)

G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=deg$Ensembl_ID,mart= mart)
empty_rows <- G_list[, 2] == ""
G_list <- G_list[!empty_rows, ]

deg2 <- merge(deg,G_list,by.x="Ensembl_ID",by.y="ensembl_gene_id", no.dups = TRUE)

deg2 <- deg2[,c(548,1:547)]
rownames(deg2) <- make.unique(as.character(deg2[, 1]))
#deg2 <- deg[,-c(1,2)]

f <- read.table("../PARG-overlap of PANoptosis and Tumor resistance genes.txt")	
deg3 <- deg2[deg2$hgnc_symbol %in% f$V1, ]

save(deg2,deg3,file = "hnsc_fpkm.rdata")

########################
####### add anno #######
########################
anno <- read.delim(file = '../TCGA-HNSC.GDC_anno.tsv',header = T)
rownames(anno) <- anno[,1]
rownames(anno) <- gsub("-", ".", rownames(anno))

#anno <- anno[,-1]
head(anno)

plotdata <- deg3[,-c(1,2)]
plotdata <- log10(plotdata)

########################
##### plot heatmap #####
########################
pheatmap::pheatmap(plotdata,scale='row',annotation_col = anno,
                   show_rownames=T,show_colnames = F,angle_col = 45)

save(plotdata,anno,file = "for_heatmap.rdata")
