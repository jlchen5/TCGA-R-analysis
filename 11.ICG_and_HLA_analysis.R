###### load packages #####

library(ggplot2)

###### load data ###### 

load("hnsc_fpkm.rdata")
ICG_and_HLA <- read.table("../ICG_and_HLA_genelist.txt",header = T)

ICG_and_HLA_exp_raw <- deg2[deg2$hgnc_symbol %in% ICG_and_HLA$gene_id , ]
ICG_and_HLA_exp_raw <- ICG_and_HLA_exp_raw[,-c(1,2)]

ICG_and_HLA_exp_anno <- t(ICG_and_HLA_exp_raw)
ICG_and_HLA_exp_anno <- merge(ICG_and_HLA_exp_anno,data.frame(subtype_df),by=0)

###### load data ######
load("ConsensusClusterPlus_results.rdata")

##### subtype info ##### 
subtype_df <- results[[2]]$consensusClass #subtype results

ICG_and_HLA_exp <- pivot_longer(ICG_and_HLA_exp_anno, 
                        cols = -c("Row.names","subtype_df"), 
                        names_to = "gene", 
                        values_to = "exp")



###### plot results
ggplot(ICG_and_HLA_exp, aes(gene,exp,fill = factor(subtype_df))) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values=c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ))+
  labs( x = "", y = "Expression") +
  stat_compare_means(method = "wilcox.test",
                     hide.ns = T,
                     size=2,
                     label = "p.format",
                     aes(group =  factor(subtype_df)))
