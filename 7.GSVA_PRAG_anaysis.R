##### load packages #####

library(GSVA)
library(msigdbr)
library("ggpubr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("clusterProfiler")
library(ggplot2)
load("cpm.rdata")
load("ConsensusClusterPlus_results.rdata")

##### !!!SET PARG Gene Sets !!! #####
hallmarks_list <- list(
  "PARG_gene_set" = c("FADD","TLR4","ITGA5","YWHAH","GSN","ITGB1",
                      "CASP6","VIM","YWHAQ","SFN","PIK3CB","DSG1",
                      "CAV1","E2F1","PIK3R2","MAPK3","PSMB10",
                      "SNAI2","PRKCD","YWHAZ"))

# First let's create a mapped data frame we can join to the gene expression values
mapped_df <- cpm
head(mapped_df)

# Map IDs from SYMBOL to ENTREZID
mapped_df2 <- data.frame(
  "entrez_id" = mapIds(
    org.Hs.eg.db,
    keys = rownames(cpm),
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
) %>%
  # Filter out NAs
  filter(!is.na(entrez_id)) %>%
  # Convert rownames to a column
  tibble::rownames_to_column("Ensembl")

save(mapped_df2,file = "symbolto_ENTREZID.rdata")


cpm_df <- merge(mapped_df,mapped_df2,by.x=0,by.y="Ensembl")

final_df <- cpm_df[,c(548:549,1:547)]

## count up how many Entrez IDs mapped to multiple Ensembl IDs
sum(duplicated(final_df$entrez_id))
rownames(final_df) <- final_df$gene_name
save(final_df,file = "cpm_with_entrez_id.rdata")


##### perform GSVA #####
gsva_results <- gsva(
  data.matrix(final_df[,-c(1:3)]),
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

score_df <- data.frame(t(gsva_results))
subtype_df <- data.frame(results[[4]]$consensusClass)
colnames(subtype_df) <- "consensusClass"
data_merged <- merge(score_df, subtype_df, by = 0)


ggplot(data_merged, aes(x = factor(data_merged$consensusClass), y = PARG_gene_set, fill = factor(data_merged$consensusClass))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  labs(title = "GSVA Scores by PARG", x = "Subtype", y = "GSVA Score") +
  theme_classic() +
  stat_compare_means(method = "wilcox.test",comparisons = list( c("1", "2"), 
                                                                c("1", "3"),
                                                                c("1", "4"),
                                                                c("2", "3"),
                                                                c("2", "4"),
                                                                c("3", "4")))




