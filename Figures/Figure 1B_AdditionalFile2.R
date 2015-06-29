#Barplot of population structure among the landraces. Figure 1B, and Additional file2

rm(list=ls())
land<-read.table("~/Documents/Structure_ready_k4.txt",header=F)

head(land)

#Add column names according to the K value used. Move external parenthesis as needed.

colnames(land)<-c("Samples","Pop1","Pop2","Pop3","Pop4","Pop5","Pop6", "Pop7")
#Select pops by majority assignment
head(land)

Pop1<-subset(land, land$Pop1 > land$Pop2 & land$Pop1 >land$Pop3 & land$Pop1 >land$Pop4 &land$Pop1 >land$Pop5 &land$Pop1 >land$Pop6 &land$Pop1 >land$Pop7)
Pop2<-subset(land, land$Pop2 > land$Pop1 & land$Pop2 >land$Pop3 & land$Pop2 >land$Pop4 & land$Pop2 >land$Pop5 & land$Pop2 >land$Pop6 & land$Pop2 >land$Pop7)
Pop3<-subset(land, land$Pop3 > land$Pop1 & land$Pop3 >land$Pop2 & land$Pop3 >land$Pop4 & land$Pop3 >land$Pop5 & land$Pop3 >land$Pop6 & land$Pop3 >land$Pop7)
Pop4 <- subset(land, land$Pop4 > land$Pop1 & land$Pop4 >land$Pop2 & land$Pop4 >land$Pop3 & land$Pop4 >land$Pop5  & land$Pop4 >land$Pop6 & land$Pop4 >land$Pop7)
Pop5 <- subset(land, land$Pop5 > land$Pop1 & land$Pop5 >land$Pop2 & land$Pop5>land$Pop3 & land$Pop5 >land$Pop4 & land$Pop5 >land$Pop6 & land$Pop5 >land$Pop7)
Pop6 <- subset(land, land$Pop6 > land$Pop1 & land$Pop6 >land$Pop2 & land$Pop6>land$Pop3 & land$Pop6 >land$Pop4 & land$Pop6 >land$Pop5 & land$Pop6 >land$Pop7)
Pop7 <- subset(land, land$Pop7 > land$Pop1 & land$Pop7 >land$Pop2 & land$Pop7>land$Pop3 & land$Pop7 >land$Pop4 & land$Pop7 >land$Pop5 & land$Pop7 >land$Pop6)

# Order each population by majority of assignment based on that population
Pop1_or<-Pop1[order(Pop1$Pop1, decreasing=T),]
Pop2_or<-Pop2[order(Pop2$Pop2, decreasing=T),]
Pop3_or<-Pop3[order(Pop3$Pop3, decreasing=T),]
Pop4_or<-Pop4[order(Pop4$Pop4, decreasing=T),]
Pop5_or<-Pop5[order(Pop5$Pop5, decreasing=T),]
Pop6_or<-Pop6[order(Pop6$Pop6, decreasing=T),]
Pop7_or<-Pop7[order(Pop7$Pop7, decreasing=T),]

# Join all the ordered populations and plot them

Pop_str<-rbind(as.data.frame(Pop1_or),as.data.frame(Pop2_or),as.data.frame(Pop3_or),as.data.frame(Pop4_or),as.data.frame(Pop5_or), as.data.frame(Pop6_or)), as.data.frame(Pop7_or)) 

pdf("~/Documents/Anita/THESIS/Data/For_STRUCTURE/Structure_land_iSelect/land_iSelect_6135_k1-k7_rep/Plot_str/iSelectLand_k5.pdf",width=10.5,height=4)
barplot(t(Pop_str[,c(2:8)]),border=0,space=0.005,col=c("purple","blue","green","brown","cyan","orange","yellow"),axisnames=F)
dev.off()  

###-------000-------- labels for when wbdc and land were analyzed together
par(mar=c(4,3,3,3))
par(xpd=T)

#Legend Figure 1A
legend(inset=-0.1,"bottom",legend=c("Central European","Coastal Mediterranean","East African","Asian"), fill=c("Blue","green","cyan","brown"),xpd=T,horiz=T)

