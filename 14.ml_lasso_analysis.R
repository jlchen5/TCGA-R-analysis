
##### load packages #####
library(glmnet)
library(survival)

##### load data #####
load("cox_raw.rdata") # df

##### lasso analysis #####
x <- as.matrix(df[,-c(1,2)])  # x为输入特征，应该是矩阵格式。否则数据框的数据会报错：Error in lognet...(串列)对象不能强制改变成'double'种类
y <- as.matrix(df$OS)
lasso <- glmnet(x = x, y = y,
                family = "binomial",
                alpha = 1,  # alpha = 1为LASSO回归，= 0为岭回归，0和1之间则为弹性网络  
                nlambda = 100)  # nlambda表示正则化路径中的个数，这个参数就可以起到一个阈值的作用，决定有多少基因的系数可以留下来。默认值为100。
print(lasso)

##### plot #####
plot(lasso, xvar = "lambda", label = TRUE)  # 系数分布图，由“log-lambda”与变量系数作图，展示根据lambda的变化情况每一个特征的系数变化，展示了Lasso回归筛选变量的动态过程
plot(lasso, xvar = "dev", label = TRUE)  # 也可以对%dev绘图
plot(lasso, xvar = "norm", label = TRUE)  # “L1范数”与变量系数作图


coef_lasso <- coef(lasso, s = 0.000555)  # 查看每个变量的回归系数。s为lasso回归的结果中的lambda值，选取相应的λ值可以得到相应数量的变量的回归系数
coef_lasso

##### cross validation #####
set.seed(1234567)  # 设置种子数，保证每次得到的交叉验证图是一致的
# 计算100个λ，画图，筛选表现最好的λ值
lasso_cv <- cv.glmnet(x = x, y = y, family = "binomial", alpha = 1, nlambda = 100)  # 交叉验证，如果结果不理想，可以重新单独运行这一行代码，或者换一下种子数
plot(lasso_cv)  # 纵坐标：部分似然偏差。上横坐标：当前的变量数。下横坐标：log(λ)
