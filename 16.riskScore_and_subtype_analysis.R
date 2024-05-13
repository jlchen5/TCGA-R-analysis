rm(list=ls())
#### riskscore and subtypes
library(ggplot2)

load("../6.无监督聚类/ConsensusClusterPlus_results.rdata")
load("riskScore_class.rdata") # riskScore_cli
subtype_df <- data.frame(results[[2]]$consensusClass) #subtype results

colnames(subtype_df) <- "subtype"

rs_sbutype_raw <- merge(subtype_df,riskScore_cli,by.x=1,by.y=0)
write.csv(rs_sbutype_raw, file = "rs_sbutype.csv")

rs_low_1 <- rs_sbutype_raw[(rs_sbutype_raw$riskScore2 == "Low" & rs_sbutype_raw$subtype == "1" ),]
rs_low_2 <- rs_sbutype_raw[(rs_sbutype_raw$riskScore2 == "Low" & rs_sbutype_raw$subtype == "2" ),]

rs_high_1 <- rs_sbutype_raw[(rs_sbutype_raw$riskScore2 == "High" & rs_sbutype_raw$subtype == "1" ),]
rs_high_2 <- rs_sbutype_raw[(rs_sbutype_raw$riskScore2 == "High" & rs_sbutype_raw$subtype == "2" ),]

subtype_count <- data_frame(
  riskclass = c("Low", "Low", "High", "High"),
  subtype = c("1", "2", "1", "2"),
  subtypecount = c(nrow(rs_low_1), nrow(rs_low_2),
                   nrow(rs_high_1), nrow(rs_high_2))
  )

write.csv(subtype_count, file = "../14.药物敏感分析/subtype_count.csv")

subtype_count %>% 
  drop_na() %>%  #去掉所有空值，避免出错
  ggplot(aes(fill=subtype, y= subtypecount, x = riskclass)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) + #纵坐标变为百分比
  labs(y = "Relative Persent %") +
  theme_bw()


