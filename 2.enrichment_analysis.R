rm(list = ls())
library(AnnotationHub)	
library(org.Hs.eg.db)   
library(clusterProfiler)
library(dplyr)
library(ggplot2)

f <- read.table("../PARG-overlap of PANoptosis and Tumor resistance genes.txt")	
f <- f[c(1)] 	#取需要的列

EG2Ensembl=toTable(org.Hs.egSYMBOL)	 #将ENTREZID和ENSEMBL对应的数据存入该变量
f=f$V1	#list转化为字符向量
geneLists=data.frame(ensembl_id=f)
results <- merge(geneLists,EG2Ensembl,by.x = "ensembl_id",by.y = "symbol",all.x=T)
id=na.omit(results$gene_id)  #提取出非NA的ENTREZID

#GO富集分析
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "ALL", pvalueCutoff = 0.01, readable= TRUE) #GO富集分析
dotplot(ego,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图

dotplot(ego, split = "ONTOLOGY", font.size = 10, #按照go中的MF/BP/CC进行换行输出
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., scale = "free") 
#scale_y_discrete(labels = function (x)str_wrap(x, width = 45))

#KEGG
ekg <- enrichKEGG(gene = id, organism = "hsa")
dotplot(ekg, font.size = 10, #按照go中的MF/BP/CC进行换行输出
        showCategory = 15)
