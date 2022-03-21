#Importing data for ANOVA
install.packages("corrplot")
setwd("D:/WOS-A2020-2023/Data_Field/19-20")
Data4AnovaY1<-read.csv("ExpData_IndMean.csv")
Data4AnovaY1[] <- lapply(Data4AnovaY1, as.numeric)
AG_dataY1<-Data4AnovaY1[,1:11]
VI_dataY1<-Data4AnovaY1[,12:29]
#Correlation amongst the investigated HVIs
corVIY1<-cor(VI_dataY1)
corrplot.mixed(round(corVIY1,2), order = 'AOE',number.cex = 0.55,tl.cex = 0.55)
#Correlation between biophysical traits and HVI
corAGVIY1<-cor(AG_dataY1,VI_dataY1)
write.csv(corAGVIY1,"corAGVIY1.csv")
library(corrplot)
corrplot(corAGVIY1, method = 'number',order = 'original',number.cex = 0.55,tl.cex = 0.55)

#PCA for data
Y1.pca<-prcomp(Data4AnovaY1,center = TRUE,scale. = TRUE)
Y1pca.rot<-as.data.frame(Y1.pca$rotation)
write.csv(Y1pca.rot, "Y1pca.csv")

desc.agy1<-describe(AG_dataY1)
desc.agy1<-desc.agy1[,c(3,4,8:13)]
#Analysis of Variance for Plant Biophysical traits
f.AGY1 = list()
p.AGY1 = list()
for(i in 2: ncol(AG_dataY1)){
  print(i)
  av<-summary(aov(DAT~AG_dataY1[,i], data = AG_dataY1))
  avv<-unclass(av)
  f.AGY1[[i-1]]<-avv[[1]][["F value"]][1]
  p.AGY1[[i-1]]<-avv[[1]][["Pr(>F)"]][1]
  i=i+1
  print(i)
}
f.AG.Y1<-as.data.frame(matrix(unlist(f.AGY1),nrow=length(f.AGY1),byrow=TRUE))
p.AG.Y1<-as.data.frame(matrix(unlist(p.AGY1),nrow=length(p.AGY1),byrow=TRUE))
f.AG.Y1<-t(f.AG.Y1)
p.AG.Y1<-t(p.AG.Y1)
fp.AG.Y1<-rbind(f.AG.Y1,p.AG.Y1)
rownames(fp.AG.Y1)<-c("f","p")
colnames(fp.AG.Y1)<-colnames(AG_dataY1[,2:21])
anov.AG.Y1<-rbind(fp.AG.Y1,AG_descdataY1)
write.csv(anov.AG.Y1,"AOV.AG.Y1.csv")

#Analysis of Variance for Hyperspectral Vegetation Indices
desc.VIy2<-describe(VI_dataY2)
desc.VIy2<-desc.VIy2[,c(3,4,8:13)]
desc.VIy2<-t(desc.VIy2)

f.VY1 = list()
p.VY1 = list()
for(i in 2: ncol(VI_dataY1)){
  print(i)
  av<-summary(aov(DAT~VI_dataY1[,i], data = VI_dataY1))
  avv<-unclass(av)
  f.VY1[[i-1]]<-avv[[1]][["F value"]][1]
  p.VY1[[i-1]]<-avv[[1]][["Pr(>F)"]][1]
  i=i+1
  print(i)
}
f.VI.Y1<-as.data.frame(matrix(unlist(f.VY1),nrow=length(f.VY1),byrow=TRUE))
p.VI.Y1<-as.data.frame(matrix(unlist(p.VY1),nrow=length(p.VY1),byrow=TRUE))
f.VI.Y1<-t(f.VI.Y1)
p.VI.Y1<-t(p.VI.Y1)
fp.VI.Y1<-rbind(f.VI.Y1,p.VI.Y1)
rownames(fp.VI.Y1)<-c("f","p")
colnames(fp.VI.Y1)<-colnames(VI_dataY1[,2:21])
anov.VI.Y1<-rbind(fp.VI.Y1,VI_descdataY1)
write.csv(anov.VI.Y1,"AOV.VI.Y1.csv")