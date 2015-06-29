#PLOT Evannon deltaK

rm(list=ls())


#Calculate DeltaK according to Evanno 2005. The format of input file (DeltaK.csv) used is:
#L(k)			L'(K)		L"(K)		DeltaK
#-1366667.67		NA		NA			NA
#-1210688.97	155978.7	67095.04	25.9795732
#-1121805.31	88883.66	26451.89	5.676267285
#-1059373.54	62431.77	19629.8		1.604109492
#-1016571.57	42801.97	7240.133333	1.389609852
#-981009.7333	35561.83667	2188.043333	0.352373022
#-947635.94		33373.79333	NA			NA
dat<-read.csv("~/Desktop/DeltaK.csv",sep=",")
pdf("~/Desktop/L.pdf",width=6,height=7)
plot(dat[,1],xlab="Clusters",ylab= c("L(k)"))
dev.off()

pdf("~/Desktop/L_primek.pdf",width=6,height=7)
plot(dat[,2],xlab="Clusters",ylab= c("L'(k)"))
dev.off()

pdf("~/Desktop/L_prime2k.pdf",width=6,height=7)
plot(dat[,3],xlab="Clusters",ylab= c("L''(k)"))
dev.off()


pdf("~/Desktop/L_Deltak.pdf",width=6,height=7)
ylab.name=expression(Delta *'K')
plot(dat[,4],xlab="Clusters",ylab=ylab.name)
lines(dat[,4], pch=22, lty=1, col="black")

dev.off()


#PLOT Clusterdness
library(plotrix)
require(plotrix)
Results_plot<-read.table("~/Clusterdness/Results_plot.txt",header=T)
pdf("~/Clusterdness/Clusterdness_plot.pdf",width=7,height=5)
plotCI(2:7,Results_plot$Average_10Runs,ui=Results_plot$UpperCI,li=Results_plot$LowerCI,ylab="Clusteredness",xlab="K clusters")
dev.off()
