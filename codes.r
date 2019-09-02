## reference from R packages
source("http://bioconductor.org/biocLite.R")
biocLite("pathifier")
biocLite("biomaRt")
biocLite("penalized")
## library the packages
library(pathifier)
library(biomaRt)
## read the dataset
GSE1456.1<-read.csv("~/Desktop/本科论文/dataset/GSE1456.csv",header=T,stringsAsFactors = F);
GSE1456<-t(GSE1456.1)
write.csv(GSE1456,"~/Desktop/本科论文/dataset/GSE1456_test.csv",seq=T,stringsAsFactors = F)
GSE1456.test<-read.csv("~/Desktop/本科论文/dataset/GSE1456_test.csv",header=T,stringsAsFactors = F)
write.table(GSE1456,"~/Desktop/本科论文/dataset/GSE1456_test1.txt")
GSE1456.anno<-read.csv("~/Desktop/本科论文/dataset/GSE1456.anno.test.csv",header=T,stringsAsFactors = F)
## 真正的数据集
a=GSE1456.1[,1]
a=c(,a)
GSE1456.test1<-read.table("~/Desktop/本科论文/dataset/GSE1456_test1.txt",header = T)
head(GSE1456.test1)
GSE1456.test2<-GSE1456.test1[,2:22284] ## 没有rownames
colnames(GSE1456.test2)<-a
## annotate the gene ID
install.packages("xml2")
library(biomaRt)
library(xml2)
listMarts() ## find out the valid Marts
ensemble=useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensemble) ## find out the dataset: 17 hsapiens_gene_ensembl
ensemble=useDataset("hsapiens_gene_ensembl",mart = ensemble)
mart=useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl") ## Testing 
c<-GSE1456.test1['ID_REF',]
attribute<-listAttributes(ensemble) ## find out the valid attributes
filters<-listFilters(ensemble) ## find out the filters
annotation=getBM(attributes=c('entrezgene','ensembl_gene_id',
                              'hgnc_symbol','affy_hg_u133_plus_2'),
                 filters="affy_hg_u133_plus_2",values=a,mart = mart)
annotation2=getBM(attributes=c('hgnc_symbol','affy_hg_u133_plus_2'),
                  filters="affy_hg_u133_plus_2",values=a,mart = mart)
## geometry mean
write.csv(annotation2,"~/Desktop/本科论文/dataset/annotation2.csv")
GSE1456.an<-cbind(annotation,)
geometry.mean <- exp(mean(log2(4)))
geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}
## transfor the annotation into excel and then combine the hgnc_symbol
## testing the average way to combine identified symbol
dat<-GSE1456.anno[19:26,]
dat<-data.frame(dat)
fun.mean<-function(x){colMeans(x[,2:ncol(dat)])}
dat.1<-t(sapply(split(dat,dat$symbol),function(x){colMeans(x[,2:ncol(dat)])}))
## 去0
newdata<-dataname[which(varname != 0),]
dataname[dataname==0]<-NA
## 去NA
new<-na.omit(dataname)
## 数值平均
GSE1456.mean<-t(sapply(split(GSE1456.anno,GSE1456.anno$symbol),function(x){colMeans(x[,2:ncol(GSE1456.anno)])}))
GSE1456.mean<-GSE1456.mean[-(1:4),] ## 去除ID为数字的gene
## 正态化检验
biocLite("nortest")
library(nortest)
lillie.test(as.matrix(GSE1456.mean)) ## 不是正态分布
## 正态化
GSE1456.mean<-data.frame(GSE1456.mean)
testdata=GSE1456.mean
colnames(testdata)<-NULL
row.names(testdata)<-NULL
testdata<-data.frame(t(GSE1456.mean))
shapiro.test(testdata)
## read the clinical dataset
GSE1456.C<-read.csv("~/Desktop/本科论文/dataset/GSE1456 Clinical dataset.csv",header=T,stringsAsFactors = F);
## Pathifier
library(pathifier)
data(KEGG)
GSE1456.mean<-data.matrix(GSE1456.mean)
Allgenes<-as.character(row.names(GSE1456.mean))
normal<-c(GSE1456.C$RELAPSE)
normal<-as.logical(normal)
pathways<-as.data.frame(read.csv("~/Desktop/本科论文/dataset/Pathways Collection.csv",stringsAsFactors = F))
pathways1<-pathways[,2:1030]
row.names(pathways1)<-pathways$X
pathways1<-as.matrix(pathways1)
install.packages("plyr")
pathways.list<-split(pathways1, row(pathways1))
pathways.list1<-split(pathways1, row(pathways1))
for(i in 1:1329)
{
  pathways1[i,][pathways1[i,]==""]=NA
  ti<-pathways1[i,][-which(is.na(pathways1[i,]))]
  n[i]<-length(ti)
}

pathways.list$`1`<-pathways.list$`1`[-which(is.na(pathways.list$`1`))]
## 把pathways.list的空格都清掉
pathways.list$`1`<-as.matrix(pathways.list$`1`[1:n[1]])
## 重复1329次上述代码
## SOM
install.packages("kohonen")
library(kohonen)
SOM.data<-t(data.matrix(PDS.data))
set.seed(1234)
pre.som<-som(scale(SOM.data),grid=somgrid(4,6,"hexagonal"))
## use hierarchical clustering to cluster the codebook vectors
groups<-40
som.hclust<-hclust(dist(scale(pre.som$distances)),method = "complete")
som.hc <- cutree(hclust(dist(pre.som$codes[[1]])), groups)
#plot
plot(pre.som, type="codes", bgcol=rainbow(groups)[som.hc])
#cluster boundaries
add.cluster.boundaries(pre.som, som.hc)
## 树图输出
plclust(som.hclust)
plot(som.hclust)
rect.hclust(som.hclust,k=40)
clustersom<-NULL
clustersom$name<-cutree(som.hclust,k=40)
u<-as.data.frame(matrix(1:13211,ncol = 1))
## 聚类后数值转换
trans.data<-matrix(t(unlist(PDS.data)),ncol=159)
trans.data<-as.data.frame(trans.data)
clustersom<-as.data.frame(clustername)
colnames(clustersom)<-c("name")
row.names(clustersom)<-c(1:13211)
trans.data1<-as.data.frame(cbind(clustersom,trans.data))
pathway.group<-split(trans.data1,trans.data1$name)
SOM.data.final<-sapply(pathway.group,function(x){colMeans(x[,2:ncol(trans.data1)])})
SOM.data.final<-as.data.frame(SOM.data.final)
row.names(SOM.data.final)<-row.names(GSE1456.test1)
colnames(SOM.data.final)<-c("G1","G2","G3","G4","G5","G6","G7","G8",
                            "G11","G10","G11","G12","G13","G14","G15","G16",
                            "G17","G18","G19","G20","G21","G22","G23","G24")
## Set up the survival models
## combine the clinical data & PDS
library(penalized)
library(survival)
clinical<-GSE1456.C[,c(3,2)]
colnames(clinical)<-c("time","event")
mydata<-cbind(clinical,SOM.data.final)
row.names(mydata)<-GSE1456.C$ID_REF
## L1-Penalized (LASSO)
prefit<-penalized(Surv(time,event),mydata[1:100,3:26],data=mydata[1:100,],lambda1 = 0.11)
coefficients(prefit) ##若是系数量级不同，可标准化。
## pretest
biocLite("globaltest")
library(globaltest)
gt(Surv(time,event), mydata[,3:42])## p value as a indicator of predictive ability
##profile the CV log likelihood
fit1 <- profL1(Surv(time,event), mydata[1:100,3:26],data = mydata[1:100,],fold=100, plot=TRUE)
## optimazing
plot(fit1$lambda, fit1$cvl, type="l")
opt1 <- optL1(Surv(time,event), mydata[1:100,3:26],data = mydata[1:100,], fold=fit1$fold)
opt1$lambda ##100  0.321319
opt1$cvl
## final model
fit<-penalized(Surv(time,event),mydata[1:100,3:26], data=mydata[1:100,],
               lambda1 = 0.24,lambda2=0)
coefficients(fit) 
## G18,重要通径簇
clustername<-as.data.frame(matrix(unlist(pre.som$unit.classif),ncol=1))
clustergene<-as.data.frame(cbind(clustername,u))
colnames(clustergene)<-c(1,2)
tPDS.data<-as.data.frame(t(PDS.data))
clustergene<-as.data.frame(cbind(clustergene,tPDS.data))
G18<-as.data.frame(clustergene[clustergene$`1`==18,])
G18<-G18[,3:161]
G18<-as.data.frame(t(G18))
mydata18<-as.data.frame(cbind(clinical,G18))
row.names(mydata18)<-GSE1456.C$ID_REF
##profile the CV log likelihood
mydata_18$time<-as.numeric(mydata_18$time)
mydata_18$event<-as.numeric(mydata_18$event)
data18<-as.data.frame(matrix(as.numeric(unlist(mydata_18)),ncol=81))
colnames(data18)<-colnames(mydata18)
row.names(data18)<-row.names(mydata18)
library(penalized)
attach(mydata18)
mydata_18<-as.matrix(mydata18)
mydata_18<-as.data.frame(mydata_18)
prefit18<-penalized(Surv(time,event),data18[1:100,3:81],data=data18[1:100,],lambda1 = 0.11)
fit18<-profL1(Surv(time,event),data18[1:100,3:81],data = data18[1:100,],fold=100, plot=TRUE) ##0.6888646
## optimazing
prefit18<-profL1(Surv(time,event), penalized=data18[1:100,3:81],data = data18[1:100,],fold=100, plot=TRUE)
plot(prefit18$lambda, prefit18$cvl, type="l")
opt18 <- optL1(Surv(time,event), data18[1:100,3:81],data = data18[1:100,], fold=100)
opt18$lambda ##100 1.374807
opt18$cvl
## final model
fit18<-penalized(Surv(time,event),data18[1:100,3:81], data=data18[1:100,],
               lambda1 = 1.374807,lambda2=0)
residuals(fit)[1:10]
fitted(fit)[1:10]
basesurv(fit)
coefficients(fit18) ##若是系数量级不同，可标准化。
loglik(fit) 
penalty(fit)
## Predict
pred<-predict(fit,mydata[140:159,3:26])
pred1<-as.data.frame(pred)
survival(pred,time=5)
plot(pred)
## Prognosis Index Source
test.Inf18<-as.matrix(data18[101:159,c(21,24,211,34,35,70,74)])
coef18<-matrix(coefficients(fit18),ncol = 1)
PIS18<-as.data.frame(test.Inf18%*%coef18)
quantile(as.vector(PIS18[,1]),probs = seq(0, 1, 0.25)) ## 确认分位数 -2.0646476 1/4   -2.1715730 1/5
data18.test<-data18[101:159,c(1,2,21,24,211,34,35,70,74)]
data18.test<-cbind(data18.test,PIS18)
data18.test<-data.frame(data18.test,rowSums(data18.test))
colnames(data18.test)<-c("time","event","V197","V241","V317","V354","V384","V1084",
                         "V1141","PIS18","HRisk")
for(i in 1:511)
{
  if(data18.test[i,10]>-2.0646476){data18.test[i,11]=c(1)}
  else{data18.test[i,11]=c(0)}
}
## KM Curve
km1456.18<-survfit(Surv(time,event)~HRisk,data = data18.test)
plot(km1456.18,col=2:3,lty=2,ylab=expression(hat(S)),
     xlab="Time",main="Observed Versus Expected Plots by Risk")
legend("bottomleft",c("LowRisk","HighRisk"),cex = 1,col=2:3,lty = 2)
## Wilcoxon log rank test
survdiff(Surv(time,event)~HRisk,data=data18.test) ##p= 0.0762 
## Predict
pred18<-predict(fit18,data18[101:159,3:81])
pred18_t<-as.data.frame(pred18@time)
S_5<-as.data.frame(survival(pred18,time=5),ncol=1)
colnames(S_5)<-c("5ys")
plot(pred18)
## 1:4比例分配
quantile(as.vector(PIS18[,1]),probs = seq(0, 1, 0.2))
for(i in 1:511)
{
  if(data18.test[i,10]>-2.1715730){data18.test[i,11]=c(1)}
  else{data18.test[i,11]=c(0)}
}
## KM Curve
km1456.18<-survfit(Surv(time,event)~HRisk,data = data18.test)
plot(km1456.18,col=2:3,lty=2,ylab=expression(hat(S)),
     xlab="Time",main="Observed Versus Expected Plots by Risk")
legend("bottomleft",c("LowRisk","HighRisk"),cex = 1,col=2:3,lty = 2)
## Wilcoxon log rank test
survdiff(Surv(time,event)~HRisk,data=data18.test) ##p= 0.0375
## ROC曲线
data18.roc<-cbind(data18.test,S_5)
install.packages("pROC")
library(pROC)
modelroc_time5<-roc(data18.roc$event,data18.roc$`5ys`)
plot(modelroc_time5, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
modelroc_Surv<-roc(data18.roc$event,data18.roc$HRisk)
plot(modelroc_Surv, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
##glmnet
##----PDS LEVEL
library(glmnet)
atime<-mydata$time[1:100]
astatus<-mydata$event[1:100]
a<-cv.glmnet(matrix(unlist(PDS.data[1:100,]),ncol=13211),
             Surv(atime,astatus),family="cox")
plot(a)
log(a$lambda.1se)
log(a$lambda.min)
b<-glmnet(matrix(unlist(PDS.data[1:100,]),ncol=13211),
          Surv(atime,astatus),family="cox",lambda = a$lambda.1se)
plot(b)
##----Gene Group
c<-cv.glmnet(matrix(unlist(mydata[1:100,3:26]),ncol=24),
             Surv(atime,astatus),family="cox")
plot(c)
c.test<-cv.glmnet(matrix(unlist(mydata[1:100,3:26]),ncol=24),
                  Surv(atime,astatus),family="cox",lambda = exp(seq(-3,-1,by=0.1)))
plot(c.test)
## 用24个Gene Group做Cox模型
fit24<-coxph(Surv(atime,astatus)~G1+G2+G3+G4+G5+G6+G7+G8+G9+G10
             +G11+G12+G13+G14+G15+G16+G17+G18+G19+G20+G21+G22+G23+G24,
             data = mydata[1:100,])
## Test
test.Inf24<-as.matrix(mydata[101:159,3:26])
coef24<-matrix(coefficients(fit24,standardize=T),ncol = 1)
PIS24<-as.data.frame(test.Inf24%*%coef24)
quantile(as.vector(PIS24[,1]),probs = seq(0, 1, 0.25)) ## 确认分位数 -8.772045 1/4   -2.1715730 1/5
data24.test<-mydata[101:159,]
data24.test<-cbind(data24.test,PIS24)
data24.test<-data.frame(data24.test,rowSums(data24.test))
colnames(data24.test)<-c("time","event","G1","G2","G3","G4","G5","G6","G7","G8","G11","G10","G11",
                         "G12","G13","G14","G15","G16","G17","G18","G19","G20","G21","G22","G23",
                         "G24","PIS24","HRisk")
for(i in 1:511)
{
  if(data24.test[i,27]>-8.772045){data24.test[i,28]=c(1)}
  else{data24.test[i,28]=c(0)}
}
## KM Curve
km1456.24<-survfit(Surv(time,event)~HRisk,data = data24.test)
plot(km1456.24,col=2:3,lty=2,ylab=expression(hat(S)),
     xlab="Time",main="Observed Versus Expected Plots by Risk")
legend("bottomleft",c("LowRisk","HighRisk"),cex = 1,col=2:3,lty = 2)
## Wilcoxon log rank test
survdiff(Surv(time,event)~HRisk,data=data24.test) ##p= 0.47 
## Predict
pred24<-predict(fit24,mydata[101:159,3:26],type=c("risk"))
pred24_t<-as.data.frame(pred24)
S_5_24<-as.data.frame(pred24,ncol=1)
colnames(S_5_24)<-c("Risk")
plot(pred24)
## ROC曲线  AUC=0.626
data24.roc<-cbind(data24.test,S_5_24)
modelroc_Surv24<-roc(data24.roc$event,data24.roc$Risk)
plot(modelroc_Surv24, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)

## 用cox来的P值来筛选最相关的
p<-c(1:24)
fitG1<-coxph(Surv(atime,astatus)~G1,data = mydata[1:100,])
fitG1
p[1]=0.907
## 重复24次上述
## penalized 下的24个Gene Group 做模型
fitPe<-penalized(Surv(time,event),mydata[1:100,3:26], data=mydata[1:100,],
               lambda1 = 0,lambda2=0)
## Prognosis Index Source
test.InfP24<-as.matrix(mydata[101:159,3:26])
coefP24<-matrix(coefficients(fitPe),ncol = 1)
PISP24<-as.data.frame(test.InfP24%*%coefP24)
quantile(as.vector(PISP24[,1]),probs = seq(0, 1, 0.25)) ## 确认分位数 -8.772045 1/4   -2.1715730 1/5
dataP24.test<-mydata[101:159,]
dataP24.test<-cbind(dataP24.test,PISP24)
dataP24.test<-data.frame(dataP24.test,rowSums(dataP24.test))
colnames(dataP24.test)<-c("time","event","G1","G2","G3","G4","G5","G6","G7","G8","G11","G10","G11",
                         "G12","G13","G14","G15","G16","G17","G18","G19","G20","G21","G22","G23",
                         "G24","PIS24","HRisk")
for(i in 1:59)
{
  if(dataP24.test[i,27]>-8.772045){dataP24.test[i,28]=c(1)}
  else{dataP24.test[i,28]=c(0)}
}
## KM Curve
km1456.P24<-survfit(Surv(time,event)~HRisk,data = dataP24.test)
plot(km1456.P24,col=2:3,lty=2,ylab=expression(hat(S)),
     xlab="Time",main="Observed Versus Expected Plots by Risk")
legend("bottomleft",c("LowRisk","HighRisk"),cex = 1,col=2:3,lty = 2)
## Wilcoxon log rank test
survdiff(Surv(time,event)~HRisk,data=dataP24.test) ##p= 0.47
## Predict
pred18<-predict(fit18,data18[101:159,3:81])
pred18_t<-as.data.frame(pred18@time)
S_5<-as.data.frame(survival(pred18,time=5),ncol=1)
colnames(S_5)<-c("5ys")
plot(pred18)

## 1:4比例分配
quantile(as.vector(PIS18[,1]),probs = seq(0, 1, 0.2))
for(i in 1:511)
{
  if(data18.test[i,10]>-2.1715730){data18.test[i,11]=c(1)}
  else{data18.test[i,11]=c(0)}
}
## KM Curve
km1456.18<-survfit(Surv(time,event)~HRisk,data = data18.test)
plot(km1456.18,col=2:3,lty=2,ylab=expression(hat(S)),
     xlab="Time",main="Observed Versus Expected Plots by Risk")
legend("bottomleft",c("LowRisk","HighRisk"),cex = 1,col=2:3,lty = 2)
## Predict
predP24<-predict(fitPe,mydata[101:159,3:26])
predP24<-as.data.frame(predP24@time)
S_P24<-as.data.frame(survival(predP24,time=5),ncol=1)
colnames(S_P24)<-c("5ys")
plot(pred18)
## Wilcoxon log rank test
survdiff(Surv(time,event)~HRisk,data=data18.test)## p= 0.0375
## ROC曲线
dataP24.roc<-cbind(dataP24.test,S_P24)
install.packages("pROC")
library(pROC)
modelroc_time5<-roc(data18.roc$event,data18.roc$`5ys`)
plot(modelroc_time5, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
modelroc_Surv<-roc(data18.roc$event,data18.roc$HRisk)
plot(modelroc_Surv, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
