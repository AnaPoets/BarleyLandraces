rm(list=ls())

#import SmartPCA output evec
DATA<-read.table("~/land_4527.evec")

DATA<-DATA[,1:803]

#Separate landrace populations. This is the output file from Clumpp_output.R
str<-read.table("~/Documents/land_assignment_k4_iSelect6135.txt",header=T)
head(str)

POP1<-subset(str,str$Populations == 1)
POP2<-subset(str,str$Populations == 2)
POP3<-subset(str,str$Populations == 3)
POP4<-subset(str,str$Populations == 4)


LAND1<-DATA[(DATA[,1] %in% POP1$Samples),]
LAND2<-DATA[(DATA[,1] %in% POP2$Samples),]
LAND3<-DATA[(DATA[,1] %in% POP3$Samples),]
LAND4<-DATA[(DATA[,1] %in% POP4$Samples),]

#pdf("~/Documents/Anita/THESIS/Data/for_smartPCA/Plots/smartPCA_1vs3.pdf",width=7,height=5)
plot(DATA[,2],DATA[,3], ylab="PC 2",xlab="PC 1",cex=0.4)
plot(DATA[,2],DATA[,3], axes=F,xlab="",ylab="",cex=0.4,col="white")
points(LAND1[,2], LAND1[,3],col="blue",cex=0.4)
points(LAND2[,2], LAND2[,3],col="red",cex=0.4)
points(LAND3[,2], LAND3[,3],col="green",cex=0.4)
points(LAND4[,2], LAND4[,3],col="cyan",cex=0.4)
abline(h=0)
abline(v=0)
dev.off()

dev.new()

#pdf("~/Documents/Anita/THESIS/Data/for_smartPCA/Plots/smartPCA_1vs2_rotated.pdf",width=7,height=5)
#plot(DATA[,3],DATA[,2], ylab="PC 2",xlab="PC 1",col="white")
plot(DATA[,3],DATA[,2], axes=F,xlab="",ylab="",cex=0.4,col="white")
points(LAND1[,3], LAND1[,2],col="blue",cex=0.4)
points(LAND2[,3], LAND2[,2],col="red",cex=0.4)
points(LAND3[,3], LAND3[,2],col="green",cex=0.4)
points(LAND4[,3], LAND4[,2],col="cyan",cex=0.4)
#abline(h=0)
#abline(v=0)
dev.off()

