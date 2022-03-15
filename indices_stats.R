#reading Indices (in tif) from directory and storing in list
files <- list.files(pattern="*.tif")
Indices <- list()
for(i in seq_along(files))
{
Indices[[i]] <- raster(files[i])
}
#saving indices list as stack
Indices_n<-stack(Indices)
#calculating statistics for indicesIn_sat
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

#normalising Indices from -1 to +1 and storing in list
norm_Indices = list()
for(i in seq_along(Indices))
{
min<-minValue(Indices[[i]])
max<-maxValue(Indices[[i]])
norm_Indices[[i]] = (Indices[[i]] - min) * (2)/(max - min)-1
}
#saving list as stack
norm_Indices_n<-stack(norm_Indices)
#calculating statistics for normalised indexes
Mean<-cellStats(norm_Indices_n,mean,na.rm=TRUE)
MD<-cellStats(norm_Indices_n,median,na.rm=TRUE)
SUM<-cellStats(norm_Indices_n,sum,na.rm=TRUE)
SD<-cellStats(norm_Indices_n,sd,na.rm=TRUE)
CV<-cellStats(norm_Indices_n,cv,na.rm=TRUE)
SK<-cellStats(norm_Indices_n,'skew',na.rm=TRUE)
stats<-data.frame("-1","1",Mean,MD,SUM,SD,CV,SK)
is.na(stats)<-sapply(stats, is.infinite)
stats[is.na(stats)]<-999
write.csv(stats,"norm_Indices_stats.csv")
setwd("../")
