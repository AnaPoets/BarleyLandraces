rm(list=ls())

Land1<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop1_contributionFromWBDC.xls",header=T)
Land2<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop2_contributionFromWBDC.xls",header=T)
Land3<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop3_contributionFromWBDC.xls",header=T)
Land4<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop4_contributionFromWBDC.xls",header=T)


SUMMARY<-matrix(ncol=5,nrow=4)

for (i in 1:4){
	landrace<-get(paste("Land",i,sep=""))
	SUMMARY[i,]<-(summary(landrace[,-c(1:4)])[4,])

}

SUMMARY_means<-sub("Mean   :","",SUMMARY)

Genome_wide<-apply(SUMMARY_means,2,as.numeric)
as.data.frame(Genome_wide)-> Genome_wide

Genome_wide_or<-Genome_wide[,c(1,2,3,4,5)]
Genome_wide_or[5,]<-apply(Genome_wide_or,2,mean)
row.names(Genome_wide_or)<-c("CentralEuropean","Asia","Coastal_Mediterranean","EastAfrica","Average")
colnames(Genome_wide_or)<-c("Northern_Mesopotamia","Southern_levant","Syrian_Desert","Northern_Levant","Central_Asia")
write.table(Genome_wide_or,"~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Genome-wide_averages.txt",quote=F,row.names=T,col.names=T,sep="\t")



#PLOT

DATA<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Genome-wide_averages.txt",header=T,row.names=1)

pdf("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Plot/Deficit_excess/Genome_wide_average.pdf",width=7,height=5)
DATA2<-DATA[,c(1,4,3,2,5)]
par(mar=c(7.5,4.5,4,7),xpd=T)
barplot(as.matrix(t(DATA2[c(1:5),])),beside=F,col=c("Red","orange","green","blue","purple"), ylab="Proportion of Ancestry",xaxt= "n",main="Genome-wide ancestry",xlab="",ylim=c(0,1))

axis(1,labels=F,tick=F)
labels=paste(c("Central European","Asian","Coastal Mediterranean","East African","Average landraces"))
text(x = c(0.65,1.9,3.1,4.3,5.5), y = par("usr")[3]-0.02, srt = 30, adj = 1,
     labels = labels, xpd = TRUE)
#legend(6.2,1,c("Central Asia","Southern Levant","Syrian Desert","Northern Mesopotamia","Northern Levant","Unassigned"),fill=c("purple","blue","green","red","orange","gray"),cex=0.6,title=expression(bold("Wild Populations")))

dev.off()





