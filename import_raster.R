#Script to import raster and perform feature extraction
setwd("D:/WOS-A2020-2023/Data_Field/20-21/HSI_20-21")
# load spatial packages
library(gplots)
library(raster)
library(rgdal)
library(RColorBrewer)
#Import the raster layer(HSI)
imgfilenames <- list.files(path="./",pattern="*.bsq",recursive=TRUE)
choice_img<-select.list(imgfilenames, preselect = NULL, multiple = FALSE,title = NULL, graphics = TRUE)
HSI_stack <- brick(choice_img)
choice_img
img_info <- file.info(choice_img)
Data_date<-as.Date(img_info$mtime)
dir_name=paste0("Data_",Data_date)
dir.create(dir_name)
#Visualising the image for finding valid outer bounds 
original_par <-par()#plot only one band
par(col.axis="white",col.lab="white",tck=0)
par(mfrow=c(1,1))
plotRGB(HSI_stack, r=70, g=50, b=30,stretch="lin", axes=TRUE, main="True color composite")
box(col="white")
mtext("Locate Area of Interest", side = 1)#locating the outer bounds on basis of user input
ext<-drawExtent()
HSI_stack_crop <- crop(HSI_stack,ext)#Crop the image stack on basis of user specifications
#Derive NDVIAparacio for masking soil from image
mask_index<-(HSI_stack_crop[[212]]-HSI_stack_crop[[112]])/(HSI_stack_crop[[212]]+HSI_stack_crop[[112]])
HSI_mask7 <- mask_index>0.7  #create binary by thresholding
HSI_canopy <-HSI_stack_crop*HSI_mask7 #extract canopy only from image
par(mfrow=c(3,1),mar=rep(1,4),oma=rep(1,4), xaxt="n", yaxt="n")
plotRGB(HSI_stack_crop, r=70, g=50, b=30,stretch="lin",axes=TRUE)
plot(HSI_mask7, col=gray.colors(2),legend=FALSE, axes=FALSE, box=FALSE)
plotRGB(HSI_canopy, r=70, g=50, b=30,stretch="lin",axes=TRUE)