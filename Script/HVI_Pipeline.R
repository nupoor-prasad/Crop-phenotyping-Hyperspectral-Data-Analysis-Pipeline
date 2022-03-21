##################################################################################
# Data Analysis Pipeline for In-Field Crop Phenotyping Using Hyperspectral Imaging
# for detection of optimum harvest time in Mentha Arvensis
#
# Nupoor Prasad, Manoj Semwal and Alok Kalra
##################################################################################
#
# set working directory
setwd("D:/Suckers_Experiment/Hyperspectral/indices")
rm(list=ls())
# load spatial packages
library(gplots)
library(raster)
library(rgdal)
library(RColorBrewer)
library(lattice)
library(corrplot)
library(psych)
library(pca3d)
library(glmnet)
#Import the hyperspectral data
imgfilenames <- list.files(path="./",pattern="*.bsq")
choice_img<-select.list(imgfilenames, preselect = NULL, multiple = FALSE,title = NULL, graphics = getOption("menu.graphics"))
HSI_stack <- brick(choice_img)
img_info <- file.info(choice_img)
Data_date<-as.Date(img_info$mtime)
dir_name=paste0("Data_",Data_date)
dir.create(dir_name)
#Visualize image for finding valid outer bounds as per user ROI 
original_par <-par()
par(col.axis="white",col.lab="white",tck=0)
plotRGB(HSI_stack, r=70, g=50, b=30,stretch="lin", axes=TRUE, main="True color composite")
box(col="white")
mtext("Locate Area of Interest", side = 1)
ext<-drawExtent()
HSI_stack_crop <- crop(HSI_stack,ext)
#Feature extraction using NDVIAparacio
mask_index<-(HSI_stack_crop[[212]]-HSI_stack_crop[[112]])/(HSI_stack_crop[[212]]+HSI_stack_crop[[112]])
HSI_mask7 <- mask_index>0.7  
HSI_canopy <-HSI_stack_crop*HSI_mask7
par(mfrow=c(3,1),mar=rep(1,4),oma=rep(1,4), xaxt="n", yaxt="n")
plotRGB(HSI_stack_crop, r=70, g=50, b=30,stretch="lin",axes=TRUE)
plot(HSI_mask7, col=gray.colors(2),legend=FALSE, axes=FALSE, box=FALSE)
plotRGB(HSI_canopy, r=70, g=50, b=30,stretch="lin",axes=TRUE)
# Calculation of hyperspectral Vegetation Indices
GNDVI <- (HSI_canopy[[145]]-HSI_canopy[[57]])/(HSI_canopy[[145]]+HSI_canopy[[57]])
HNDVI <- (HSI_canopy[[179]]-HSI_canopy[[109]])/(HSI_canopy[[179]]+HSI_canopy[[109]])
RENDVI <- (HSI_canopy[[145]]-HSI_canopy[[125]])/(HSI_canopy[[145]]+HSI_canopy[[125]])
NDVIAparacio <- (HSI_canopy[[212]]-HSI_canopy[[114]])/(HSI_canopy[[212]]+HSI_canopy[[114]])
NDVIHaboudane <- (HSI_canopy[[167]]-HSI_canopy[[110]])/(HSI_canopy[[167]]+HSI_canopy[[110]])
NDVIZarcoTejada <- (HSI_canopy[[156]]-HSI_canopy[[113]])/(HSI_canopy[[156]]+HSI_canopy[[113]])
NDVIrouse<-(HSI_canopy[[194]]-HSI_canopy[[110]])/(HSI_canopy[[194]]+HSI_canopy[[110]])
PSNDa <- (HSI_canopy[[167]]-HSI_canopy[[114]])/(HSI_canopy[[167]]+HSI_canopy[[114]])
PSNDb <- (HSI_canopy[[167]]-HSI_canopy[[94]])/(HSI_canopy[[167]]+HSI_canopy[[94]])
CIRe <- (HSI_canopy[[167]]/HSI_canopy[[134]])-1
MCARI1 <- 1.2*(2.5*(HSI_canopy[[167]]-HSI_canopy[[110]])-1.3*(HSI_canopy[[167]]-HSI_canopy[[57]]))
PSRI <- (HSI_canopy[[113]] - HSI_canopy[[35]]) / (HSI_canopy[[145]])
PRI <- (HSI_canopy[[48]] - HSI_canopy[[66]]) / (HSI_canopy[[48]] + HSI_canopy[[66]])
Datt<-HSI_canopy[[110]]/(HSI_canopy[[126]]*HSI_canopy[[57]])
BGBO<-(HSI_canopy[[33]]/HSI_canopy[[59]])-(HSI_canopy[[33]]/HSI_canopy[[94]])
ARI<-((1/HSI_canopy[[57]])-(1/HSI_canopy[[123]]))*HSI_canopy[[167]]
CRIr<-((1/HSI_canopy[[39]])-(1/HSI_canopy[[123]]))*HSI_canopy[[167]]
CRIg<-((1/HSI_canopy[[39]])-(1/HSI_canopy[[57]]))*HSI_canopy[[167]]
MSR <- (HSI_canopy[[167]]-HSI_canopy[[11]])/(HSI_canopy[[114]]-HSI_canopy[[11]])
SIPI <- (HSI_canopy[[167]] - HSI_canopy[[11]]) / (HSI_canopy[[167]] - HSI_canopy[[114]])
Indices<-stack(GNDVI,HNDVI,RENDVI,NDVIAparacio,NDVIHaboudane,NDVIZarcoTejada,NDVIrouse,PSNDa,PSNDb,CIRe,MCARI1,PSRI,PRI,Datt,BGBO,ARI,CRIr,CRIg,SIPI,MSR)
names(Indices)<-c("GNDVI","HNDVI","RENDVI","NDVIAparacio","NDVIHaboudane","NDVIZarcoTejada","NDVIrouse","PSNDa","PSNDb","CIRe","MCARI1","PSRI","PRI","Datt","BGBO","ARI","CRIr","CRIg","SIPI","MSR")
#Visualise rasters and histograms of indices
par(mfrow=c(7,6),mai = c(0.3, 0.2, 0.1, 0.1))
for( i in 1:20 ) 
{
plot( Indices[[i]],axes=FALSE,box=FALSE)
hist(Indices[[i]],axes=FALSE,box=FALSE,xlab=NULL,ylab=NULL)
title(names(Indices[[i]]))
}
Sys.sleep(5)
setwd(dir_name)
writeRaster(Indices, filename=names(Indices), bylayer=TRUE,format="GTiff",overwrite=TRUE)
setwd("../")
#calculating statistics for indices
files <- list.files(pattern="*.tif")
Indices <- list()
for(i in seq_along(files))
{
Indices[[i]] <- raster(files[i])
}
Indices_n<-stack(Indices)
Min<-cellStats(Indices_n,'min',na.rm=TRUE)
Max<-cellStats(Indices_n,'max',na.rm=TRUE)
Mean<-cellStats(Indices_n,mean,na.rm=TRUE)
MD<-cellStats(Indices_n,median,na.rm=TRUE)
SUM<-cellStats(Indices_n,sum,na.rm=TRUE)
SD<-cellStats(Indices_n,sd,na.rm=TRUE)
CV<-cellStats(Indices_n,cv,na.rm=TRUE)
SK<-cellStats(Indices_n,'skew',na.rm=TRUE)
stats1<-data.frame(Min,Max,Mean,MD,SUM,SD,CV,SK)
is.na(stats1)<-sapply(stats1, is.infinite)
stats[is.na(stats1)]<-999
write.csv(stats1,"Indices_stats.csv")
#Assimilating stats of all indices over all dates
files <- list.files("./", recursive = T, pattern = "norm_Indices_stats.csv")
In_stats <- list()
for(i in seq_along(files))
{
In_stats[[i]] <- read.csv(files[i])
}
#calculating basic image statistics for all indexes through the acquisition dates
MN <- lapply(In_stats, function(x) x$Mean)
MD <- lapply(In_stats, function(x) x$MD)
SD <- lapply(In_stats, function(x) x$SD)
CV <- lapply(In_stats, function(x) x$CV)
SK <- lapply(In_stats, function(x) x$SK)
Mean_stats <- do.call(cbind, MN)
Median_stats <- do.call(cbind, MD)
Stdv_stats <- do.call(cbind, SD)
CVar_stats <- do.call(cbind, CV)
Skew_stats <- do.call(cbind, SK)
L <- list(Mean_stats,Median_stats,Stdv_stats,CVar_stats,Skew_stats)
names(L)<-c("Mean_stats","Median_stats","Stdv_stats","CVar_stats","Skew_stats")
In_names<-In_stats[[1]]$X
L <- lapply(L, function(x){colnames(x) <- c("156DAS","165DAS","176DAS","187DAS","197DAS")
 							x
				  }
		)
L <- lapply(L, function(x){rownames(x) <- In_names
							x
				  }
		)
lapply(1:length(L), function(i) write.csv(L[[i]],file = paste0(names(L[i]), ".csv")))

#Importing data (ground measurements plus image indices) for statistical analysis
setwd("D:/WOS-A2020-2023/Data_Field/19-20")
Data4AnovaY1<-read.csv("ExpData_IndMean_new.csv")
Data4AnovaY1[] <- lapply(Data4AnovaY1, as.numeric)
setwd("D:/WOS-A2020-2023/Data_Field/20-21")
Data4AnovaY2<-read.csv("Exp_Mean_new.csv")
Data4AnovaY2[] <- lapply(Data4AnovaY2, as.numeric)
#Calculating descriptive statistics for biophysical traits and hyperspectral indices
#Year1
desc_dataY1<-describe(Data4AnovaY1)
desc_dataY1<-desc_dataY1[,c(3,4)]
desc_dataY1$cv<-(desc_dataY1[,2]/desc_dataY1[,1])*100
desc_dataY1<-t(desc_dataY1)
#Year2
desc_dataY2<-describe(Data4AnovaY2)
desc_dataY2<-desc_dataY2[,c(3,4)]
desc_dataY2$cv<-(desc_dataY2[,2]/desc_dataY2[,1])*100
desc_dataY2<-t(desc_dataY2)
#Anova for biophysical traits and hyperspectral indices
#year1
f.Y1 = list()
p.Y1 = list()
for(i in 2: ncol(Data4AnovaY1)){
aov<-summary(aov(DAT~Data4AnovaY1[,i], data = Data4AnovaY1))
aov<-unclass(aov)
f.Y1[[i-1]]<-aov[[1]][["F value"]][1]
p.Y1[[i-1]]<-aov[[1]][["Pr(>F)"]][1]
i=i+1
}
f.Y1<-as.data.frame(matrix(unlist(f.Y1),nrow=length(f.Y1),byrow=TRUE))
p.Y1<-as.data.frame(matrix(unlist(p.Y1),nrow=length(p.Y1),byrow=TRUE))
f.Y1<-t(f.Y1)
p.Y1<-t(p.Y1)
fp.Y1<-rbind(f.Y1,p.Y1)
rownames(fp.Y1)<-c("f","p")
colnames(fp.Y1)<-colnames(Data4AnovaY1[,-1])
#year2
f.Y2 = list()
p.Y2 = list()
for(i in 2: ncol(Data4AnovaY2)){
  aov<-summary(aov(DAT~Data4AnovaY2[,i], data = Data4AnovaY2))
  aov<-unclass(aov)
  f.Y2[[i-1]]<-aov[[1]][["F value"]][1]
  p.Y2[[i-1]]<-aov[[1]][["Pr(>F)"]][1]
  i=i+1
}
f.Y2<-as.data.frame(matrix(unlist(f.Y2),nrow=length(f.Y2),byrow=TRUE))
p.Y2<-as.data.frame(matrix(unlist(p.Y2),nrow=length(p.Y2),byrow=TRUE))
f.Y2<-t(f.Y2)
p.Y2<-t(p.Y2)
fp.Y2<-rbind(f.Y2,p.Y2)
rownames(fp.Y2)<-c("f","p")
colnames(fp.Y2)<-colnames(Data4AnovaY2[,-1])
#combining results of descriptive stats and ANOVA for evaluation
#Year1
ANOVA_Y1<-rbind(desc_dataY1[,-1],fp.Y1)
setwd("D:/WOS-A2020-2023/Data_Field/19-20")
write.csv(ANOVA_Y1, 'ANOVA_Y1.csv')
#Exploring Correlation amongst the investigated Biophysical Traits
corY1<-cor(Data4AnovaY1)
corrplot.mixed(round(corY1[2:11,2:11],2), order = 'original',number.cex = 0.55,tl.cex = 0.55)
write.csv(corY1[2:11,2:11],"corAGY1.csv")
#Correlation between biophysical traits and HVI
corAGVIY1<-cor(Data4AnovaY1[,2:11],Data4AnovaY1[,12:29])
write.csv(corAGVIY1,"corAGVIY1.csv")
corrplot(corAGVIY1, method = 'number',order = 'original',number.cex = 0.55,tl.cex = 0.55)
#Year2
ANOVA_Y2<-rbind(desc_dataY2[,-1],fp.Y2)
setwd("D:/WOS-A2020-2023/Data_Field/20-21")
write.csv(ANOVA_Y2, 'ANOVA_Y2.csv')
#Exploring Correlation amongst the investigated Biophysical Traits
corY2<-cor(Data4AnovaY2)
corrplot.mixed(round(corY2[2:11,2:11],2), order = 'original',number.cex = 0.55,tl.cex = 0.55)
write.csv(corY2[2:11,2:11],"corAGY2.csv")
#Correlation between biophysical traits and HVI
corAGVIY2<-cor(Data4AnovaY2[,2:11],Data4AnovaY2[,12:29])
write.csv(corAGVIY2,"corAGVIY2.csv")
corrplot(corAGVIY2, method = 'number',order = 'original',number.cex = 0.55,tl.cex = 0.55)
#PCA for biophysical traits and hyperspectral indices in Y1
#Year1
Data4PCAY1<-Data4AnovaY1[,-1]
norm_Data4PCAY1 = list()
for(i in 1: ncol(Data4PCAY1)){
min<-min(Data4PCAY1[[i]])
max<-max(Data4PCAY1[[i]])
norm_Data4PCAY1[[i]] = (Data4PCAY1[[i]] - min) * (2)/(max - min)-1
i=i+1
}
norm_Data4PCAY1<-as.data.frame(norm_Data4PCAY1)
names(norm_Data4PCAY1)<-names(Data4PCAY1)
Y1.pca<-prcomp(norm_Data4PCAY1,center = TRUE,scale. = TRUE)
Y1pca.rot<-as.data.frame(Y1.pca$rotation)
PCA_Y1<-rbind(Y1pca.rot,summary(Y1.pca)$importance)
sig_pca_rot<-subset(Y1pca.rot[,1:3], PC1>0.25|PC2>0.25|PC3>0.25|PC1<(-0.25)|PC2<(-0.25)|PC3<(-0.25))#eliminating insignificant values
Sig_PCA_Y1<-rbind(sig_pca_rot,summary(Y1.pca)$importance[,1:3])
write.csv(Sig_PCA_Y1, "SIGPCA_Y1.csv")
pca3d(Y1.pca, components = 1:3,biplot=TRUE, biplot.vars=4,show.plane=FALSE)
snapshotPCA3d(file="pca3dY1.png")
#Year2
Data4PCAY2<-Data4AnovaY2[,-1]
norm_Data4PCAY2 = list()
for(i in 1: ncol(Data4PCAY2)){
  min<-min(Data4PCAY2[[i]])
  max<-max(Data4PCAY2[[i]])
  norm_Data4PCAY2[[i]] = (Data4PCAY2[[i]] - min) * (2)/(max - min)-1
  i=i+1
}
norm_Data4PCAY2<-as.data.frame(norm_Data4PCAY2)
names(norm_Data4PCAY2)<-names(Data4PCAY2)
Y2.pca<-prcomp(norm_Data4PCAY2,center = TRUE,scale. = TRUE)
Y2pca.rot<-as.data.frame(Y2.pca$rotation)
PCA_Y1<-rbind(Y2pca.rot,summary(Y2.pca)$importance)
sig_pca_rot2<-subset(Y2pca.rot[,1:3], PC1>0.25|PC2>0.25|PC3>0.25|PC1<(-0.25)|PC2<(-0.25)|PC3<(-0.25))#eliminating insignificant values
Sig_PCA_Y2<-rbind(sig_pca_rot2,summary(Y2.pca)$importance[,1:3])
write.csv(Sig_PCA_Y1, "SIGPCA_Y1.csv")
pca3d(Y2.pca, components = 1:3,biplot=TRUE, biplot.vars=4,show.plane=FALSE)
snapshotPCA3d(file="pca3dY2.png")
#Lasso for biophysical traits and hyperspectral indices
Data4LASSOY1Y2<-rbind(Data4AnovaY1,Data4AnovaY2)
order_Data4LASSOY1Y2<-Data4LASSOY1Y2[ order(Data4LASSOY1Y2$DAT) , ]
order_Data4LASSOY1Y2<-order_Data4LASSOY1Y2[,c(4:14,16:22,24:29)]
norm_order_Data4LASSOY1Y2 = list()
for(i in 1: ncol(order_Data4LASSOY1Y2)){
min<-min(order_Data4LASSOY1Y2[[i]])
max<-max(order_Data4LASSOY1Y2[[i]])
norm_order_Data4LASSOY1Y2[[i]] = (order_Data4LASSOY1Y2[[i]] - min) * (2)/(max - min)-1
i=i+1
}
norm_order_Data4LASSOY1Y2 <-as.data.frame(norm_order_Data4LASSOY1Y2)
names(norm_order_Data4LASSOY1Y2)<-names(order_Data4LASSOY1Y2)
rsq=list()
Std_Beta=list()
rmse=list()
x <- data.matrix(norm_order_Data4LASSOY1Y2[,c(9:24)])
AllSD = apply(norm_order_Data4LASSOY1Y2[,c(9:24)],2,sd)
for(i in 1:8){
j=i*100
set.seed(j)
y<-norm_order_Data4LASSOY1Y2[,i]
cv_model <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
y_predicted <- predict(best_model, s = best_lambda, newx = x)
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)
rsq[[i]] <- 1 - sse/sst
Std_Beta[[i]] = coefficients(best_model)[-1]*AllSD
mse = function(x,y) {mean((x-y)^2)}
mse(y, y_predicted)
rmse[[i]]<-mse(y, y_predicted)
i=i+1
}
rsq<-as.data.frame(rsq)
names(rsq)<-colnames(norm_order_Data4LASSOY1Y2[,1:8])
Std_Beta<-as.data.frame(Std_Beta)
names(Std_Beta)<-colnames(norm_order_Data4LASSOY1Y2[,1:8])
rmse<-as.data.frame(rmse)
names(rmse)<-colnames(norm_order_Data4LASSOY1Y2[,1:8])
Mod_perf<-rbind(rsq,rmse)
rownames(Mod_perf)<-c("rsq","rmse")
plot(test[,1]~factor(rownames(test)),type="p",col="red",pch=16,ylab="Value",xlab="Biophysical Trait")
points(test[,2],col="red",pch=16)
title(main = "Model performance for biophysical traits", cex.main = 1)
legend(0.3, 1, legend=c("rmse", "rsq"),col=c("red", "black"), lty=1:1, cex=0.8)
write.csv(Mod_perf, "model_performance.csv")