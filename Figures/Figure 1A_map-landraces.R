rm(list=ls())
#Input str2 is the Q matrix from STRUCTURE
str2<-read.table("~/Documents/output_names_str_k4.txt",header=F)
head(str2)
#colnames(str2)<-c("Samples","Pop1","Pop2","Pop3","Pop4","Pop5","Pop6", "Pop7")
colnames(str2)<-c("Samples","Pop1","Pop2","Pop3","Pop4")
head(str2)
#Call land long locations for all the samples
latlong_land<-read.table("~/Documents/Genetic_map.txt",header=T)
head(latlong_land)
#Add header names
colnames(latlong_land)<-c("Samples", "HABIT", "SPIKEROW","COUNTRY","LATITUDE","LONGITUDE", "ELEV.m.","Noise_data", "LATITUDE_original", "LONGITUDE_original")
head(latlong_land)
#Order by sample name
latlong_land[order(latlong_land$Samples),]->latlong_land_order
head(latlong_land_order)
dim(latlong_land_order)


###latlong coordinates with a little bit of noise so the dots don't overlap that much
lonlat_corrected<-read.table("~/Documents/803_land_CORRECT_latlong.txt",header=T)
head(lonlat_corrected)
latlong_land_order<-lonlat_corrected[order(lonlat_corrected $Samples),]
#########

#Pick up the samples that are in the Structure
latlong_land_order[(latlong_land_order$Samples %in% str2[,1]),]->latlong_land_str
str2[order(str2[,1]),]->str_order
setequal(latlong_land_str[,1],str_order[,1])
head(str_order)

#colnames(str_order)<-c("Samples","Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7")
head(str_order)
#Combine structure results with latlong coordinates
#cbind(as.data.frame(str_order), as.data.frame(latlong_land_str[,c(5,6)]))->dat
cbind(as.data.frame(str_order), as.data.frame(latlong_land_str[,c(3,4)]))->dat
dim(dat)
head(dat)



##PLOT PIE CHARTS
library(maps)
library(mapdata)
library(mapplots)
library(scales)

#pdf("~/Desktop/Land_pops_prop.pdf",width=7,height=5)
map("worldHires",xlim = c(-9,145),ylim=c(5,70),col=gray(c(0.6)),fill=F)

#map("worldHires",xlim = c(-9,145),ylim=c(5,70),fill=T,col="lemonchiffon3")

for (i in 1:(dim(dat)[1])){
add.pie(z=as.numeric(dat[i,2:5]), x=dat[i,7], y=dat[i,6], col=c("blue","red","green","cyan"),radius=sqrt(0.6), labels="")
}
#dev.off()
