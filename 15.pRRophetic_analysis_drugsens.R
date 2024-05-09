####### load packages ####### 
#install.packages("oncoPredict",binary = T)
library(oncoPredict)
load("hnsc_fpkm.rdata")

####### load data ####### 

GDSC2_Expr <- readRDS(file="./DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS(file = "./DataFiles/Training Data/GDSC2_Res.rds")

testExpr <- deg2[,-c(1,2)]
testExpr <- log2(testExpr+1)

####### 计算所有药物敏感性 #######

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = as.matrix(testExpr),#需要matrix
              batchCorrect = 'eb',  
              #IC50是对数转换的，所以表达矩阵也用对数转换过的
              powerTransformPhenotype = F,
              minNumSamples = 20,
              printOutput = T,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = "homogenizeData"
)

## 结果就一个文件，就是每个样本对每一个药物的IC50值，读取进来看看：
res <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
dim(res)
## [1] 414 546
res[1:4,1:4]

## 有了这个结果，我们就可以取出感兴趣的药物，可视化IC50值在不同组间的差异,这里随便取前10个药物：
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

res %>% 
  select(1:11) %>% 
  bind_cols(sample_group = sample_group) %>% 
  pivot_longer(2:11,names_to = "drugs",values_to = "ic50") %>% 
  ggplot(., aes(sample_group,ic50))+
  geom_boxplot(aes(fill=sample_group))+
  scale_fill_jama()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(vars(drugs),scales = "free_y",nrow = 2)+
  stat_compare_means()


