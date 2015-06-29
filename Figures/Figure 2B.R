#This script is used to plot results from supportmix
# 0. Make a table with the proportion of ancestry contributed by all wild in each landrace. Use the mean values to set the Zero line for Delta ancestry
# 1. Delta ancestry in different panels
# 2. Delta ancestry all in one panel
# 3. Chromosomal painting proportions

rm(list=ls())

POP1<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop1_contributionFromWBDC.xls",header=T)
POP2<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop2_contributionFromWBDC.xls",header=T)
POP3<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop3_contributionFromWBDC.xls",header=T)
POP4<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop4_contributionFromWBDC.xls",header=T)

GENOME_WIDE_PROPORTIONS<-matrix(ncol=5,nrow=5)
GENOME_WIDE_PROPORTIONS [1,]<-c(mean(POP1[,5]),mean(POP1[,6]),mean(POP1[,7]),mean(POP1[,8]),mean(POP1[,9]))
GENOME_WIDE_PROPORTIONS [2,]<-c(mean(POP2[,5]),mean(POP2[,6]),mean(POP2[,7]),mean(POP2[,8]),mean(POP2[,9]))
GENOME_WIDE_PROPORTIONS [3,]<-c(mean(POP3[,5]),mean(POP3[,6]),mean(POP3[,7]),mean(POP3[,8]),mean(POP3[,9]))
GENOME_WIDE_PROPORTIONS [4,]<-c(mean(POP4[,5]),mean(POP4[,6]),mean(POP4[,7]),mean(POP4[,8]),mean(POP4[,9]))
AverageGenomeWide<-(apply(GENOME_WIDE_PROPORTIONS[1:4,],2,mean))
GENOME_WIDE_PROPORTIONS[5,]<-AverageGenomeWide

row.names(GENOME_WIDE_PROPORTIONS)<-c("Central European","Asia","Coastal Mediterranean","East Africa","Average")
colnames(GENOME_WIDE_PROPORTIONS)<-c("Northern_Mesopotamia","Southern_Levant","Syrian_Desert","Northern_Levant","Central_Asia")
write.table(GENOME_WIDE_PROPORTIONS,"~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Genome_wide_proportions.xls",quote=F,row.names=T,col.names=T,sep="\t")

#Call EACH LANDRACE !!population with its contribution from the wild

fred <-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/Pop3_contributionFromWBDC.xls",header=T)
colnames(fred)<-c("Interval","Chro","index","Cumulative","Northern_Mesopotamia","Southern_Levant","Syrian_Desert","Northern_Levant","Central_Asia")
summary(fred)

# Using only 5 wild populations summary(fred)[4,]), for all pops and take the average
#yel<-0.1959475 
#oran<-0.110338 
#red<-0.086175 
#green<-0.05905325 
#blue<-0.456245
#pink<-0.092287 

oran<-GENOME_WIDE_PROPORTIONS[5,4]
red<-GENOME_WIDE_PROPORTIONS[5,1]
green<-GENOME_WIDE_PROPORTIONS[5,3]
blue<-GENOME_WIDE_PROPORTIONS[5,2]
pink<-GENOME_WIDE_PROPORTIONS[5,5]

#(fred$Unassigned - yel )/yel->fred$Unassigned
(fred$Northern_Levant - oran)/oran ->fred$Northern_Levant
(fred$Northern_Mesopotamia - red )/red->fred$Northern_Mesopotamia
#(fred$Caspian_Sea - darkg )/darkg-> fred$Caspian_Sea
(fred$Syrian_Desert -green)/green -> fred$Syrian_Desert
(fred$Southern_Levant - blue)/blue ->fred$Southern_Levant
(fred$Central_Asia -pink)/pink -> fred$Central_Asia


## 1. Plot all pops wild in one pannel
plot(c(0,1112.71),c(min(for_plot[,6:10])+0.02,15), col="white", ylab="", xlab="",xaxt='n')
minimum<-(round(min(for_plot[,6:10])))-1

#for east africa and asia
#plot(c(0,1112.71),c(min(for_plot[,6:11])+0.02,20), col="white",ylab="",yaxt="n", xlab="",xaxt='n')
#axis(2, pos=c(minimum), at=c(minimum,0,5,10,20))
#axis(1, pos=c(minimum), at=c(0,200,400,600,800,1000,1113))
#segments(0,20,1113,20)
#segments(1113, c(minimum),1113,20)


points(for_plot$Positions,for_plot$Northern_Levant,cex=0.1,col="Orange")
points(for_plot$Positions,for_plot$Northern_Mesopotamia,cex=0.1,col="red")
#points(for_plot$Positions,for_plot$Caspian_Sea,cex=0.1,col="dark green")
points(for_plot$Positions,for_plot$Syrian_Desert,cex=0.1,col="green")
points(for_plot$Positions,for_plot$Southern_Levant,cex=0.1,col="blue")
points(for_plot$Positions,for_plot$Central_Asia,cex=0.1,col="pink")
#points(for_plot$Positions,for_plot$Unassigned,cex=0.1,col="yellow")


##lty 1,2,3,4,5 to have different type of liine per population
#lines(for_plot$Positions,for_plot$Unassigned,col="yellow")
lines(for_plot$Positions,for_plot$Central_Asia,col="purple",lty=1,lwd=2)
lines(for_plot$Positions,for_plot$Northern_Levant, col ='orange',lty=1,lwd=2)
lines(for_plot$Positions,for_plot$Northern_Mesopotamia,col="red",lty=1,lwd=2)
#lines(for_plot$Positions,for_plot$Caspian_Sea,col="dark green",lty=1,lwd=2)
lines(for_plot$Positions,for_plot$Syrian_Desert,col="green",lty=1,lwd=2)
lines(for_plot$Positions,for_plot$Southern_Levant,col="blue",lty=1,lwd=2)

#lines to devide chromosomes
abline(v=c(146.9,325.23,489.65,620.45,804.38,943.77),col=c(grey(0.50)),xpd=FALSE,lty=5,lwd=0.5)

#significant deltaAncestry as the extrem values of the excess/deficit of all wild pops the whole genome
head(for_plot)
odd_th<-seq(1, nrow(for_plot), by=2)
data<-for_plot[c(odd_th),]
threshold<-quantile(as.matrix(data[,c(6:10)]),0.95)
## 2.475145 all populations together

abline(h= threshold,lty=3)


#abline(v=c(fred$Cumulative[2:dim(fred)[1]]), col=c(grey(0.97)), xpd=F)

# xaxis linkage group name at the median position of the cM distance per chromosome
#axis(side=1,at=c(55.30,218.3,392.1,545.5, 703.2,863.6, 1022.0), labels=c("1H","2H","3H","4H","5H","6H","7H"),)
#Genome-wide
abline(h=0,col="black",lty=5,lwd=1)

#Legend outside margin
#legend(734,15.6, c("Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia"),col=(c("orange","red","green","blue","pink")),cex=1,lty=1,lwd=2,bg="white",xpd=T,title=expression(bold("Wild Populations"))) 

legend("topleft", c("Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia", "98th percentile","Genome-wide average ancestry"),col=(c("orange","red","green","blue","purple","black","black")),cex=0.97,lty=c(rep(1,5),3,5),lwd=c(rep(2,5),1,1),bg="white",title=expression(bold("Wild Populations")))

dev.off()
#Special SNPs with annotations
#axis(side=3, at=c(236.81,287.52,528.09,550.37,559.98,620.55,663.25,746.26,864.09), labels=c("Vpe2d/ Vrs1","Hbc14/ Pap6", "PhyA1/ Pap20","Agl4/ S85","","Hina/Gsp1","TaSAG3","VRN-1","Cry2"),las=3,cex.axis=0.45)
axis(side=3, at =c(69.88,622.99,746.26,756.55), labels=c("Nec1","GSP/ Hina", "HKO1/ VRN-H1", "Dhn9"),las=3,cex.axis=0.45)

dev.off()



## 2 . PLOT chromosomal painting proportions 

rbind(Pop_summary[,5],Pop_summary[,6],Pop_summary[,7],Pop_summary[,8],Pop_summary[,9],Pop_summary[,10])->VALUES
dim(VALUES)
#We want bars with different widths representing the cM in each chromosomal segment
Pop_summary$Interval ->width

#pdf(file="~/Dropbox/ANITA_301213/THESIS/Data/for_SupporMix/Supportmix_barley1846_ACTG_TO_AB/Wild_k5/75_windows/plots/Pop1_proportions.pdf", width=8, height=5)

barplot(VALUES,width=0.3,beside=FALSE,col=c("red","blue","green","orange","purple","yellow"),space=0)

