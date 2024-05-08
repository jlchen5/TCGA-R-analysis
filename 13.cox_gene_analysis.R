rm(list=ls())
###### load packages ######
library("survival")
library("survminer")
library(tidyr)


##### load data #####
load("surival_and_subtype.rdata") # surival data: km_raw
load("hnsc_fpkm.rdata") # gene expression data: deg2 and deg3(n=20)

km <- km_raw[,]
exp_data_raw <- deg3[,-2]
exp_data <- t(exp_data_raw[,-1])

##### merge them #####
cox_raw <- merge(km_raw,exp_data,by.x="sample",by.y=0)
rownames(cox_raw) <- cox_raw[,1]
df <- cox_raw[,-c(1,2,4)]
save(df,file = "./cox_raw.rdata")

##### cox analysis #####
#设置p值的阈值
pfilter <- 0.05  
#新建空白数据框
uniresult <- data.frame()  
#使用for循环对输入数据中的100个基因依次进行单因素COX分析
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(df[,3:ncol(df)])){   
  unicox <- coxph(Surv(time = OS.time, event = OS) ~ df[,i], data = df)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  if(pvalue<pfilter){ 
    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  }
}   
#保存单因素COX回归分析结果
write.csv(uniresult,file = "single_factor_COX_results.csv",row.names = F)


###### plot ######

tducs <- read.csv("single_factor_COX_results.csv",header = T)
rownames(tducs) <- tducs$"gene"
tducs <- tducs[,-1]
gene <- rownames(tducs)
hr <- tducs$"HR"
hrLow  <- tducs$"L95CI"
hrHigh <- tducs$"H95CI"
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.05, "<0.05", sprintf("%s", tducs$pvalue))

pdf(file="UniCoxSurForestPlot.pdf", width = 6,height = 3)
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))

plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)

dev.off()

