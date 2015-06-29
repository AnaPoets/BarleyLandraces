#Make supplemental table 5. Proportions of ancestry and great circle distance from the closest wild population
rm(list=ls())

c_asia<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_C_Asia.txt",header=T,row.names=1)

n_levant<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_N_levant.txt",header=T,row.names=1)

n_mesopotamia<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_N_mesopotamia.txt",header=T,row.names=1)

s_desert<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_S_desert.txt",header=T,row.names=1)

s_levant<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_S_levant.txt",header=T,row.names=1)

dim(c_asia)
#Grab the ancestry and distance from each population. eg. if Asia, grab Central_Asia and Distance columns from the input file c_asia

Proportions_dis<-cbind(as.data.frame(n_levant[,6]),as.data.frame(n_mesopotamia[,3]),as.data.frame(s_desert[,5]),as.data.frame(s_levant[,4]),as.data.frame(c_asia[,7]),as.data.frame(n_levant[,8]),as.data.frame(n_mesopotamia[,8]),as.data.frame(s_desert[,8]),as.data.frame(s_levant[,8]),as.data.frame(c_asia[,8]))

row.names(Proportions_dis)<-row.names(c_asia)
colnames(Proportions_dis)<-c("Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia","Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia")


head(Proportions_dis)

#Round to two decimals

Proportions_round<-round(Proportions_dis,digits=2)
head(Proportions_round)

#Separate by landrace populations according to the genetic assignment obtained from Structure
str<-read.table("~/Documents/land_assignment_k4_iSelect6135.txt",header=T)
head(str)

land1<-subset(str,str$Populations ==1)
land2<-subset(str,str$Populations ==2)
land3<-subset(str,str$Populations ==3)
land4<-subset(str,str$Populations ==4)

Prop_land1<-Proportions_round[(row.names(Proportions_round) %in% land1$Samples),]
Prop_land2<-Proportions_round[(row.names(Proportions_round) %in% land2$Samples),]
Prop_land3<-Proportions_round[(row.names(Proportions_round) %in% land3$Samples),]
Prop_land4<-Proportions_round[(row.names(Proportions_round) %in% land4$Samples),]

EMPTY_ROW<-as.data.frame(t(rep("",(dim(Proportions_round)[2]))))
colnames(EMPTY_ROW)<-c("Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia","Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia")


Ind_prop_distance<-rbind(as.data.frame(Prop_land1),as.data.frame(EMPTY_ROW),as.data.frame(Prop_land2),as.data.frame(EMPTY_ROW),as.data.frame(Prop_land3),as.data.frame(EMPTY_ROW),as.data.frame(Prop_land4))

#write.table(Ind_prop_distance,"~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/TableS5_prop_distance.xls",quote=F,row.names=T,col.names=T,sep="\t")

#PLOT BOXPLOT
#Grab only percentage of assignment, no distance
Pop1<-Prop_land1[,1:5]
Pop2<-Prop_land2[,1:5]
Pop3<-Prop_land3[,1:5]
Pop4<-Prop_land4[,1:5]

#Plot barplots of the distribution of proportion of ancestry at individuals from each landrace population

boxplot(Pop1[,1],Pop2[,1],Pop3[,1],Pop4[,1],Pop1[,2],Pop2[,2],Pop3[,2],Pop4[,2],Pop1[,3],Pop2[,3],Pop3[,3],Pop4[,3],Pop1[,4],Pop2[,4],Pop3[,4],Pop4[,4],Pop1[,5],Pop2[,5],Pop3[,5],Pop4[,5],col=c("blue","Red","Green","Cyan","blue","Red","Green","Cyan","blue","Red","Green","Cyan","blue","Red","Green","Cyan","blue","Red","Green","Cyan"))

#Using ggplot

Pop1_new<-cbind(rep("Pop1",(dim(Pop1)[1])),Pop1)
colnames(Pop1_new)<-c("Population","Northern Levant", "Northern Mesopotamia", "Syrian Desert", "Southern Levant","Central Asia")

Pop2_new<-cbind(rep("Pop2",(dim(Pop2)[1])),Pop2)
colnames(Pop2_new)<-c("Population", "Northern Levant", "Northern Mesopotamia", "Syrian Desert", "Southern Levant","Central Asia")

Pop3_new<-cbind(rep("Pop3",(dim(Pop3)[1])),Pop3)
colnames(Pop3_new)<-c("Population", "Northern Levant", "Northern Mesopotamia", "Syrian Desert", "Southern Levant","Central Asia")


Pop4_new<-cbind(rep("Pop4",(dim(Pop4)[1])),Pop4)
colnames(Pop4_new)<-c("Population", "Northern Levant", "Northern Mesopotamia", "Syrian Desert", "Southern Levant","Central Asia")


DATA_pops<-rbind(Pop1_new, Pop2_new, Pop3_new, Pop4_new)
head(DATA_pops)

require(reshape2)

DATA_transformed<-melt(DATA_pops, id.var= "Population")


library(ggplot2)

pdf("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Plot/Distribution_proportions_ancestry.pdf",width=11 , height=5)
bp<-ggplot(data=DATA_transformed, aes(x=variable, y =value, fill= Population)) + geom_boxplot()

bp+ scale_fill_manual(name= "Landrace Populations", values=c("blue","red","green","cyan"), labels=c("Central European","Asian","Coastal Mediterranean","East African")) + labs(x="Wild populations",y="Proportion Ancestry")


dev.off()


