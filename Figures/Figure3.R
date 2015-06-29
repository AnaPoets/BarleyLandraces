rm(list=ls())
#call the population assignement according to Structure k=4.
POPULATIONS<-c("C_asia","N_levant","N_mesopotamia","S_desert","S_levant")

for (p in 1:length(POPULATIONS)){

prop_latlong<-read.table(paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_",POPULATIONS[p],".txt",sep=""),header=T)

#Change column names so they match with POPULAITONS
colnames(prop_latlong)<-c("SAMPLES","LATITUDE_original", "LONGITUDE_original" ,"N_mesopotamia", "S_levant" ,"S_desert", "N_levant" ,"C_asia"  ,"Distance")
assi<-read.table("~/Dropbox/ANITA_301213/THESIS/Data/For_STRUCTURE/Structure_w_out_learning_samples_landraces/land_assignment_k4_iSelect6135.txt",header=T)

# Select samples names for each population

pop1<-subset(assi,assi$Populations == "1")
pop2<-subset(assi,assi$Populations == "2")
pop3<-subset(assi,assi$Populations == "3")
pop4<-subset(assi,assi$Populations == "4")


# Get only individuals growing in the native range of wild barley
prop_latlong <-subset(prop_latlong, prop_latlong$LATITUDE_original >= 29.72 & prop_latlong$LATITUDE_original <=46.77 & prop_latlong$LONGITUDE_original >=20.9 & prop_latlong $LONGITUDE_original <=73.66)

# Separate resutls from supportMix for each population

prop_latlong[(prop_latlong[,1] %in% pop1$Samples),]->supmix1
prop_latlong[(prop_latlong[,1] %in% pop2$Samples),]->supmix2
prop_latlong[(prop_latlong[,1] %in% pop3$Samples),]->supmix3
prop_latlong[(prop_latlong[,1] %in% pop4$Samples),]->supmix4

head(prop_latlong)


##PLOT just the summary of the proportion ancestry distribution from each wild population
ave_dis1<-median(supmix1$Distance)
ave_dis2<-median(supmix2$Distance)
ave_dis3<-median(supmix3$Distance)
ave_dis4<-median(supmix4$Distance)
ave_distance<-rbind(ave_dis1,ave_dis2,ave_dis3,ave_dis4)

#PLOT each of the landrace populations

pdf(paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Plot/Ancestry_vs_distance/",POPULATIONS[p],"_ave_wRange.pdf",sep=""),width=7,height=5)

HEADER<-c("Central Asia","Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant")
YAXIS<-c(0.5,0.5,0.5,0.5,1)
boxplot(supmix1[,POPULATIONS[p]],xlim=c(0,2050),ylim=c(0, YAXIS[p]),at=c(ave_dis1),boxwex = 150,col="blue", main=(paste(HEADER[p]," wild ancestry", sep="")), ylab="Proportion Ancestry",xlab="Great circle distance (Km)")
boxplot(supmix2[,POPULATIONS[p]],add=T,at=ave_dis2,boxwex=150,col="red")
boxplot(supmix3[,POPULATIONS[p]],add=T,at=ave_dis3,boxwex=150,col="green")

axis(1, labels =T)

if (p == 2){
legend("topright",col=c("green","blue","red","cyan"),fill=c("Green","Blue","Red","Cyan"),legend=c("Coastal Mediterranean","Central European","Asian"),cex=0.75,bg="white",title=expression(bold("Landraces Populations")))
	}

dev.off() 

##Plot with correlation frequency vs distance. CHANGE COLUMN OF POPULATION ANALYZED!!
print(cor(prop_latlong$Distance,prop_latlong[,POPULATIONS[p]]))
}