#Format CLUMPP output to a more informative file, adding sample names and removing extra columns
#This code takes Clumpp output format and fix it to plot the Structure results for each K value

rm(list=ls())

for (k in 1:7){
#Take one of the Structure outputs to extract the sample names.
# Example of Structure output format
# CIho00497  7 0.0002 0.0003 0.9993 0.0001
# CIho01458  7 0.6556 0.1895 0.0867 0.0682
# CIho01461  7 0.6581 0.1542 0.1495 0.0383
# CIho01468  7 0.5973 0.2087 0.1015 0.0925
# CIho01604  7 0.0001 0.9998 0.0000 0.0001
# CIho02266  7 0.7899 0.0064 0.0573 0.1464

Structure_out<-read.table("~/Documents/Landrace_str_K4_R1_q",header=F)


#Call CLUMPP output file
#Example of CLUMPP output file:
#1  1 (1)  8  : 0.0001 0.0002 0.9993 0.0003
#2  2 (1)  8  : 0.0677 0.6557 0.0874 0.1892
#3  3 (1)  8  : 0.0390 0.6577 0.1489 0.1544
#4  4 (1)  8  : 0.0920 0.5974 0.1015 0.2090
#5  5 (1)  8  : 0.0001 0.0001 0.0001 0.9997
#6  6 (1)  8  : 0.1459 0.7901 0.0581 0.0059

Clumpp_out<-read.table(paste("~/CLUMPP/output_k",K,"/Land_eu_af_asi_iSelect6135_k",k,".outfile",sep=""),header=F)

Structure_clumpp<-cbind(Structure_out[,1],Clumpp_out[,-c(1:5)])


write.table(Structure_clumpp,paste("~/Documents/Structure_ready_k",k,".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

}


#For k4 match the samples with the lat long coordinates for future analysis

latlong<-read.table("~/Documents/github/BarleyLandraces/Datasets/Samples_latlong.txt",header=T)
#Select the 803 samples used
latlong_803<-latlong[(latlong$SAMPLES %in% Structure_clumpp[,1]),]
#Order by sample name
latlong_803_or<-latlong_803[order(latlong_803$SAMPLES),]

if(identical(as.character(latlong_803_or$SAMPLES),as.character(Structure_clumpp[,1])) == FALSE) stop (print("Sample names in different order"))

land <-cbind(Structure_clumpp [,1],latlong_803_or[,c(5,6)], Structure_clumpp[,-c(1)])
colnames(land)<-c("Samples","LATITUDE","LONGITUDE","Pop1","Pop2","Pop3","Pop4")
#Select pops by majority assignment

Pop1<-subset(land, land$Pop1 > land$Pop2 & land$Pop1 >land$Pop3 & land$Pop1 >land$Pop4 )
Pop2<-subset(land, land$Pop2 > land$Pop1 & land$Pop2 >land$Pop3 & land$Pop2 >land$Pop4)
Pop3<-subset(land, land$Pop3 > land$Pop1 & land$Pop3 >land$Pop2 & land$Pop3 >land$Pop4)
Pop4 <-subset(land, land$Pop4 > land$Pop1 & land$Pop4 >land$Pop2 & land$Pop4 >land$Pop3)

Pop_order<-rbind(Pop1,Pop2,Pop3,Pop4)

Land_assignment_k4_iSelect6152<-cbind(Pop_order[,c(1:3)],c(rep(1,(dim(Pop1)[1])),rep(2,(dim(Pop2)[1])),rep(3,(dim(Pop3)[1])),rep(4,(dim(Pop4)[1]))))

colnames(Land_assignment_k4_iSelect6152)<-c("Samples","LATITUDE","LONGITUDE","Populations")

write.table(Land_assignment_k4_iSelect6152,"~/Documents/Land_assignment_k4_iSelect6152.txt",quote=F,row.names=F,col.names=T,sep="\t")