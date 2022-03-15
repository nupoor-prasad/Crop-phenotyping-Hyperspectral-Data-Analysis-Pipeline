# calculate Indices
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
#Visualise Rasters and histogram
par(mfrow=c(7,6),mai = c(0.3, 0.2, 0.2, 0.1))
for( i in 1:20 ) 
{
plot( Indices[[i]],axes=FALSE,box=FALSE)
hist(Indices[[i]],axes=FALSE,box=FALSE,xlab=NULL,ylab=NULL)
title(names(Indices[[i]]))
}
Sys.sleep(5)
dev.copy(jpeg,paste0("indices_",Data_date,".jpg"),width=500, height=600, res=70)
dev.off()
setwd(dir_name)
# Export raster
writeRaster(Indices, filename=names(Indices), bylayer=TRUE,format="GTiff",overwrite=TRUE)
setwd("../")


