##### load packages ##### 
library(survminer) 
library(survival) 

###### load data ######
load("ConsensusClusterPlus_results.rdata")

##### subtype info ##### 
subtype_df <- results[[4]]$consensusClass #subtype results
subtype_df <- data.frame(subtype_df)
table(subtype_df)

##### surivial info ##### 
surivial_df <- read.delim("../TCGA-HNSC.survival.tsv",header = T)
surivial_df$sample <- gsub(pattern = "-",replacement = ".",surivial_df$sample)
head(surivial_df)

##### merge them ##### 
km_raw <- merge(subtype_df,surivial_df,by.x=0,by.y="sample")
colnames(km_raw)=c("sample","subtype","OS","X_PATIENT","OS.time")
head(km_raw)



km_data <- survfit(Surv(time = km_raw$OS.time, event = km_raw$OS)~subtype,
                                            data=km_raw)
##### plot Survival cruve
ggsurvplot(km_data ,
           pval = T,
           legend.title = "subtype_df",
           legend="right",
           #font.main = c(16, "bold", "darkblue"),
           #censor.shape=26,
           #palette = c("blue","red"),
           censor=F,
           xlab="overall survival",
           title="survival analysis")

##### plot heatmap for cluster 
load("for_heatmap.rdata")

### add subtype to anno
merge_anno <- merge(subtype_df,anno,by=0)
rownames(merge_anno) <- merge_anno$Row.names
merge_anno$subtype_df <- factor(merge_anno$subtype_df)
merge_anno <- merge_anno[,c(2,7,9)]


pheatmap::pheatmap(plotdata,scale='row',
                   annotation_col = merge_anno,
                   show_rownames=T,show_colnames = F,angle_col = 45)
