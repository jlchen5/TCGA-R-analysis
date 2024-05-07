###### TME Anlaysis with ESTIMATE and CIBERSORT ######

###### ESTIMATE first ###### 

##### load package #####
#install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
library(estimate)
library(ggplot2)
library("ggpubr")

###### load data  #####
load("hnsc_fpkm.rdata") #deg2(all fpkm) and deg3(PARG fpkm)
load("ConsensusClusterPlus_results.rdata")
subtype_df <- data.frame(results[[4]]$consensusClass)
colnames(subtype_df) <- "subtype"


esti_raw <- deg2[,-c(1,2)]
write.table(esti_raw,file="./exp.txt",sep = "\t",quote = F)

filterCommonGenes(input.f = "./exp.txt",
                  output.f = "esti_filt.gct",
                  id = "GeneSymbol")

estimateScore(input.ds="./esti_filt.gct",
              output.ds="./esti_score.gct",
              platform = "illumina")

scores=read.table("esti_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores

TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])
scores_new <- data.frame(scores,TumorPurity)
scores_merge <- merge(scores_new,subtype_df,by=0)


##
ggplot(scores_merge, aes(x = factor(scores_merge$subtype), y = scores_merge$ImmuneScore, fill = factor(scores_merge$subtype))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  labs(title = " ", x = "Subtype", y = "ImmuneScore Scores") +
  theme_classic()+
  stat_compare_means(method = "wilcox.test",comparisons = list( c("1", "2"), 
                                                                c("1", "3"),
                                                                c("1", "4"),
                                                                c("2", "3"),
                                                                c("2", "4"),
                                                                c("3", "4")))


###### then CIBERSORT ###### 
##### load packages #####
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(preprocessCore)
library(CIBERSORT)
library(tidyr)

##### set data
load("deg_res.rdata")
LM22.file <- read.delim("./LM22.txt")
rownames(LM22.file) <- LM22.file[,1]
LM22.file <- data.matrix(LM22.file[,-1])
TCGA_exp.file <- data.matrix(deg2[,-c(1,2)])


#####
TCGA_TME.results <- cibersort(LM22.file ,TCGA_exp.file, perm = 50, QN = F)  



TME_data <- data.frame(TCGA_TME.results[,1:22])
TME_data_anno <- merge(TME_data,subtype_df,by=0)

TME_New <- pivot_longer(TME_data_anno, 
                        cols = -c(Row.names, subtype), 
                        names_to = "cell_type", 
                        values_to = "value")



ggplot(TME_New, aes(cell_type,value,fill = factor(subtype))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  labs(title = "TME Cell composition", x = "Celltype", y = "GSVA Score") +
  theme_classic() 




