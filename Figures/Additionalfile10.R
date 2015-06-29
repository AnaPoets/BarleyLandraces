##PLOT
#After runing Step3_SupportMix_calculate_proportions_wild_testing_1896.R
#for with unassigned
dir.create(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/plot",sep=""))

pdf(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/plot","/Run",D,"_proportions_w_unassigned.pdf",sep=""),width=9,height=5)
par(mar=c(4,4,4,8),xpd=FALSE)
barplot(as.matrix(t(SUMMARY_w_unassigned_plot)),col=c("red","blue","green","orange","purple","gray"),ylim=c(0,1),xaxt="n"); 
axis(1,at=c(.7,1.88,3.1,4.3,5.5),lab=c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"),cex.axis=.65)
legend(6.1,1, c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"),fill=(c("red","blue","green","orange","purple")),cex=0.7,expression(bold(title="Wild Populations")),xpd=T)

dev.off()

#for with unassigned

pdf(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/plot","/Run",D,"_proportions_all_assigned.pdf",sep=""),width=9,height=5)
par(mar=c(4,4,4,8),xpd=FALSE)
barplot(as.matrix(t(SUMMARY_all_assigned_plot)),col=c("red","blue","green","orange","purple"),ylim=c(0,1),xaxt="n"); 
axis(1,at=c(.7,1.88,3.1,4.3,5.5),lab=c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"),cex.axis=.65)
legend(6.1,1, c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"),fill=(c("red","blue","green","orange","purple")),cex=0.7,expression(bold(title="Wild Populations")),xpd=T)

dev.off()

}

######=================================NOW GET TOGETHER ALL THE RESULTS IN A NICE SUMMARY!!!=================================================
##SUMMARY BOXPLOT ACCROSS 50 RESPS


#Colect the results for each wild testing population, accross all the 50 reps

#Create directories to place summaries
SUMMARY_PLOTS<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,"Plots_summary",sep="")
dir.create(SUMMARY_PLOTS)

SUMMARY_PLOTS_W_UNASSIGNED<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,"Plots_summary/","WithUnassigned",sep="")
dir.create(SUMMARY_PLOTS_W_UNASSIGNED)


SUMMARY_PLOTS_ALL_ASSIGNED <-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,"Plots_summary/","All_assigned",sep="")
dir.create(SUMMARY_PLOTS_ALL_ASSIGNED)

#Choose input file 1) Genome_wide_prop_summary_w_unassigned or 2) Genome_wide_prop_summary_all_assigned
# 1) Genome_wide_prop_summary_w_unassigned  
RESULTS_N_mesop_ALL<-matrix(ncol=6,nrow=FINALIZE)
RESULTS_S_Levant<-matrix(ncol=6,nrow= FINALIZE)
RESULTS_SyrianDesert<-matrix(ncol=6,nrow= FINALIZE)
RESULTS_N_Levant<-matrix(ncol=6,nrow= FINALIZE)
RESULTS_C_Asia<-matrix(ncol=6,nrow= FINALIZE)

for (i in  INITIATE:FINALIZE){


infile<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,"20_wild_testing_rep",i,"/Analysis/Pop_summary/","Genome_wide_prop_summary_w_unassigned.txt",sep="")
myInfile<-read.table(infile,header=T)
RESULTS_N_mesop_ALL[i,1:6]<-as.matrix(myInfile[1,])
RESULTS_S_Levant[i,1:6]<-as.matrix(myInfile[2,])
RESULTS_SyrianDesert[i,1:6]<-as.matrix(myInfile[3,])
RESULTS_N_Levant[i,1:6]<-as.matrix(myInfile[4,])
RESULTS_C_Asia[i,1:6]<-as.matrix(myInfile[5,])
}

LOCATION<-SUMMARY_PLOTS_W_UNASSIGNED

pdf(paste(SUMMARY_PLOTS_W_UNASSIGNED,"/N_Mesopotamia.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_N_mesop_ALL[,1:6],main="Northern mesopotamia (testing samples)",xaxt= "n", xlab=" ",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia","Unassigned"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_W_UNASSIGNED,"/S_Levant.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_S_Levant[,1:6],main="Southern Levant (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia","Unassigned"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_W_UNASSIGNED,"/S_Desert.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_SyrianDesert[,1:6],main="Syrian Desert (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia","Unassigned"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_W_UNASSIGNED,"/N_levant.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_N_Levant[,1:6],main="Northern Levant (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia","Unassigned"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_W_UNASSIGNED,"/C_Asia.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_C_Asia[,1:6],main="Central Asia (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia","Unassigned"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

MEAN_SUMMARY<-matrix(nrow=5,ncol=6)
row.names(MEAN_SUMMARY)<-c("N_Mesopotamia","S_Levant","Syrian Desert","Nothern_Levant","Central_Asia")
colnames(MEAN_SUMMARY)<-c("N_Mesopotamia","S_Levant","Syrian Desert","Nothern_Levant","Central_Asia","Unassigned")
MEAN_SUMMARY[1,]<-apply(RESULTS_N_mesop_ALL,2,mean)
MEAN_SUMMARY[2,]<-apply(RESULTS_S_Levant,2,mean)
MEAN_SUMMARY[3,]<-apply(RESULTS_SyrianDesert,2,mean)
MEAN_SUMMARY[4,]<-apply(RESULTS_N_Levant,2,mean)
MEAN_SUMMARY[5,]<-apply(RESULTS_C_Asia,2,mean)

#Write table with summary of proportion of ancestry.
Output_table<-paste(SUMMARY_PLOTS_W_UNASSIGNED,"/Summary_proportion_Ancestry_50_run_w_unassigned.xls",sep="")
write.table(MEAN_SUMMARY, Output_table,quote=F,row.names=T,col.names=T,sep="\t")


2) Genome_wide_prop_summary_all_assigned

RESULTS_N_mesop_ALL2<-matrix(ncol=5,nrow=FINALIZE)
RESULTS_S_Levant2<-matrix(ncol=5,nrow= FINALIZE)
RESULTS_SyrianDesert2<-matrix(ncol=5,nrow= FINALIZE)
RESULTS_N_Levant2<-matrix(ncol=5,nrow= FINALIZE)
RESULTS_C_Asia2<-matrix(ncol=5,nrow= FINALIZE)


for (i in  INITIATE:FINALIZE){


infile2<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,"20_wild_testing_rep",i,"/Analysis/Pop_summary/","Genome_wide_prop_summary_all_assigned.txt",sep="")
myInfile2<-read.table(infile2,header=T)
RESULTS_N_mesop_ALL2[i,1:5]<-as.matrix(myInfile2[1,])
RESULTS_S_Levant2[i,1:5]<-as.matrix(myInfile2[2,])
RESULTS_SyrianDesert2[i,1:5]<-as.matrix(myInfile2[3,])
RESULTS_N_Levant2[i,1:5]<-as.matrix(myInfile2[4,])
RESULTS_C_Asia2[i,1:5]<-as.matrix(myInfile2[5,])
}



pdf(paste(SUMMARY_PLOTS_ALL_ASSIGNED,"/N_Mesopotamia.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_N_mesop_ALL2[,1:5],main="Northern mesopotamia (testing samples)",xaxt= "n", xlab=" ",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_ALL_ASSIGNED,"/S_Levant.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_S_Levant2[,1:5],main="Southern Levant (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_ALL_ASSIGNED,"/S_Desert.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_SyrianDesert2[,1:5],main="Syrian Desert (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_ALL_ASSIGNED,"/N_levant.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_N_Levant2[,1:5],main="Northern Levant (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste(SUMMARY_PLOTS_ALL_ASSIGNED,"/C_Asia.pdf",sep=""),width=7,height=6)
par(mar=c(11,4,4,4))
boxplot(RESULTS_C_Asia2[,1:5],main="Central Asia (testing samples)",xaxt= "n", xlab="",ylab="Proportion of assignment",ylim=c(0:1))
axis(1,labels=F)
labels=paste(c("Northern Mesopotamia","Southern Levant","Syrian Desert","Northern Levant","Central Asia"))
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

MEAN_SUMMARY2<-matrix(nrow=5,ncol=5)
row.names(MEAN_SUMMARY2)<-c("N_Mesopotamia","S_Levant","Syrian Desert","Nothern_Levant","Central_Asia")
colnames(MEAN_SUMMARY2)<-c("N_Mesopotamia","S_Levant","Syrian Desert","Nothern_Levant","Central_Asia")
MEAN_SUMMARY2[1,]<-apply(RESULTS_N_mesop_ALL2,2,mean)
MEAN_SUMMARY2[2,]<-apply(RESULTS_S_Levant2,2,mean)
MEAN_SUMMARY2[3,]<-apply(RESULTS_SyrianDesert2,2,mean)
MEAN_SUMMARY2[4,]<-apply(RESULTS_N_Levant2,2,mean)
MEAN_SUMMARY2[5,]<-apply(RESULTS_C_Asia2,2,mean)

#Write table with summary of proportion of ancestry.
Output_table<-paste(SUMMARY_PLOTS_ALL_ASSIGNED,"/Summary_proportion_Ancestry_50_run_all_assigned.xls",sep="")
write.table(MEAN_SUMMARY2, Output_table,quote=F,row.names=T,col.names=T,sep="\t")