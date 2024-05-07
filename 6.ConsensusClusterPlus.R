##### Load Packages #####
library(ConsensusClusterPlus)

##### Load data #####
load("cpm.rdata")
TumorMat <- cpm[,-547]

results <- ConsensusClusterPlus(as.matrix(TumorMat), 
                                maxK = 6, reps = 200, 
                                pItem = 0.8,
                                pFeature = 0.8, 
                                clusterAlg = "hc", 
                                distance = "pearson", 
                                innerLinkage="complete",
                                title = "ConsensusCluster_results", 
                                plot = "pdf")
save(results, file= "ConsensusClusterPlus_results.rdata")
# 计算cluster-consensus and item-consensus

icl = calcICL(results,title="icl_results",plot="pdf")
#icl是两元素的列表，分别是cluster-consensus和item-consensus测度
icl[["clusterConsensus"]]
results[[4]]$consensusClass #subtype results




