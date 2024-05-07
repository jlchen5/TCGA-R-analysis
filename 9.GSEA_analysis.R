######## GSEA_analysis #####

###### load packages ######
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

###### load data ######
load("deg_res.rdata")
load("symbolto_ENTREZID.rdata") # which is mapped_df2
load("cpm_with_entrez_id.rdata") # which is final_df

###### sort the genelist ###### 
head(res)
res_anno_EntrezID <- merge(res,mapped_df2,by.x=0,by.y="Ensembl")
# 按照 log2FoldChange 列降序排序整个数据框
res_anno_EntrezID <- res_anno_EntrezID[order(res_anno_EntrezID$log2FoldChange, decreasing = TRUE), ]

# 创建 geneList
gene_list <- setNames(res_anno_EntrezID$log2FoldChange, res_anno_EntrezID$entrez_id)

# 确保 geneList 是降序排列的
gene_list <- sort(gene_list, decreasing = TRUE)


###### 执行 GSEA ######
gsea_result <- gseKEGG(geneList = gene_list, 
                       organism = "hsa", 
                       #nPerm = 1000, 
                       #minGSSize = 10, 
                       #maxGSSize = 500, 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       verbose = FALSE)

####### 查看结果 ######
head(summary(gsea_result))


gseaplot2(gsea_result, geneSetID = 1:5)

# 使用 dotplot 可视化基因集的 P 值和富集得分
dotplot(gsea_result, showCategory = 20) + ggtitle("Dotplot of Gene Sets")

# 使用 cnetplot 显示富集的基因集及其相关基因的网络图
cnetplot(gsea_result, foldChange = gene_list) + ggtitle("Network Plot of Gene Sets and Genes")

