#Plot each landrace population : frequency of each wild private allele vs genetic map
#After Steps 1 and 2 from ~/github/BarleyLandraces/Analysis/Shared_private have run
rm(list=ls())

#Land1 = Central European, Land2= Asia, Land3=Coastal Mediterranean , Land4= East Africa
CentralEuropean<-read.table("~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/LAND_1_private_freq.txt",header=T)
Asia<-read.table("~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/LAND_2_private_freq.txt",header=T)
CoastalMediterranean<-read.table("~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/LAND_3_private_freq.txt",header=T)
EastAfrica<-read.table("~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/LAND_4_private_freq.txt",header=T)

#Separate wild populations
#Select which landrace population to work with
land_pops<-c("CentralEuropean","Asia","CoastalMediterranean","EastAfrica")
for (i in c(1,3)){	
POPULATION<- get(paste(land_pops[4]))

n_mesopotamia<-subset(POPULATION,POPULATION$Population == 'N_Mesopotamia')
s_levant<-subset(POPULATION,POPULATION$Population == 'S_Levant')
s_desert<-subset(POPULATION,POPULATION$Population == 'S_Desert')
n_levant<-subset(POPULATION,POPULATION$Population == 'N_Levant')
c_asia<-subset(POPULATION,POPULATION$Population == 'C_Asia')

#PLOT

OUTPUT<-paste("~/Documents/plots/",land_pops[i],".pdf",sep="")
pdf(OUTPUT,width=7,height=5) #asia with labels width=8
#par(mar=c(4,4,1,12),xpd=T)
#plot(POPULATION$Cumulative,POPULATION$frequency_private,xaxt='n',ylab='Frequency',xlab="Linkage Group",cex=0.5,ylim=c(0,1),col="white")

#for Asia and East Africa with no Y-AXIS
plot(POPULATION$Cumulative,POPULATION$frequency_private,xaxt='n',yaxt='n',xlab="Linkage Group",ylab="",cex=0.5,ylim=c(0,1),col="white")

points(s_levant$Cumulative, s_levant$frequency_private,col='Blue',cex=0.8)
points(s_desert$Cumulative, s_desert$frequency_private,col='green',cex=0.8)
points(n_levant$Cumulative, n_levant$frequency_private,col='orange',cex=0.8)
points(c_asia$Cumulative, c_asia$frequency_private,col='magenta',cex=0.8)
points(n_mesopotamia$Cumulative, n_mesopotamia$frequency_private,col='red',cex=0.8)

axis(side=1,at=c(55.30,218.3,392.1,545.5, 705.2,863.6, 1022.0), labels=c("1H","2H","3H","4H","5H","6H","7H"))
#lines to devide chromosomes
abline(v=c(146.9,325.23,503.65,620.45,813.38,935.77),col=c(grey(0.50)),xpd=FALSE,lty=5,lwd=0.5)

#legend(1190,1.04,col=c("red","blue","green","orange","magenta"),pch=c(1,1,1,1),legend=c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"),cex=0.97,bg="white",title=expression(bold("Wild Populations")))

dev.off()

}