rm(list=ls())

LIST_SEGMENTS<-read.table("~/Documents/PLINK/Sep_30SNP/output/Analysis/List_segment_w_matches.txt",header=F)

#Get segments info and genetic map for 30SNPs exact division, created with Step2_Plink_IBD_1896_input_map_ped_30SNP.R
SEGMENT_information<-read.table("~/Documents/PLINK/Sep_30SNP/Segments_info.txt")


genmap_30SNPs<-read.table("~/Documents/PLINK/Sep_30SNP/map_cm_30snpALL.map")
genmap_iSelect<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)


#1.Loop through all the output files from PLINK after trimming from significants. Count how many wild and landrace share IBD by each segment.
#2.Add a Segment identifiyer and put all segments together for plotting

CountWildLand<-matrix(ncol=2,nrow=(dim(LIST_SEGMENTS)[1]))
row.names(CountWildLand)<-LIST_SEGMENTS[,1]
colnames(CountWildLand)<-c("Wild","Land")

ave_wild<-as.data.frame((mean(CountWildLand[,1]))/277) #277 is the total number of wild samples after removing the Caspian Sea population
ave_land<-as.data.frame((mean(CountWildLand[,2]))/803) # 803 is the total number of landraces



for (i in 1:(dim(LIST_SEGMENTS)[1])){
	
	SEGMENT_name<-paste(LIST_SEGMENTS[i,1],"_perfectMatch.txt",sep="")
	SEGMENT<-read.table(paste("~/Documents/PLINK/Sep_30SNP/output/Analysis/",SEGMENT_name,sep=""))
	
	#find number of unique wild and landrace in this segment, either from column 2 or 4
	SAMPLES<-c(as.character(SEGMENT[,2]), as.character(SEGMENT[,4]))
	SAMPLES_unique<-sort(unique(SAMPLES[order(SAMPLES)]))
	CountWildLand[i,1]<-length(grep("^WBDC", SAMPLES_unique))
	CountWildLand[i,2]<-length(SAMPLES_unique)-length(grep("^WBDC", SAMPLES_unique))
	
	#Add a column identifier to each segment file. Save the new file in a new variable
	IDENTIFIER<-rep(LIST_SEGMENTS[i,1],(dim(SEGMENT)[1]))
	SEGMENT_ID<-cbind(as.data.frame(IDENTIFIER),as.data.frame(SEGMENT))
	SEGMENT_ID_IMPORTANT<-as.data.frame(SEGMENT_ID[,c(1,3,5)])
	assign(paste(LIST_SEGMENTS[i,1],"_important",sep=""),SEGMENT_ID_IMPORTANT)
	}

#Combine all the segments with IDs

SegmentsID_names<-grep(x=ls(pos=1),pattern="_important",value=T)

All_important_segments<-do.call(what=rbind,args=mget(SegmentsID_names))

#Average cM distance between start and end at each chromosomal segment

Segments_contributed_atLeast_ONECE<-SEGMENT_information[(SEGMENT_information[,1] %in% All_important_segments$IDENTIFIER),]

RANGE<-(Segments_contributed_atLeast_ONECE[,11]-Segments_contributed_atLeast_ONECE[,5])

RESULTS_SUMMARY<-matrix(ncol=2,nrow=6)
RESULTS_SUMMARY[1,1]<-c("average_wild")
RESULTS_SUMMARY[1,2]<-ave_wild[1,1]
RESULTS_SUMMARY[2,1]<-c("average_land")
RESULTS_SUMMARY[2,2]<-ave_land[1,1]
RESULTS_SUMMARY[3,1]<-c("Min cM")
RESULTS_SUMMARY[3,2]<-c(min(RANGE))
RESULTS_SUMMARY[4,1]<-c("Max cM")
RESULTS_SUMMARY[4,2]<-c(max(RANGE))
RESULTS_SUMMARY[5,1]<-c("Mean cM")
RESULTS_SUMMARY[5,2]<-c(mean(RANGE))
RESULTS_SUMMARY[6,1]<-c("No. of segments that have ~100% match")
RESULTS_SUMMARY[6,2]<-c(dim(Segments_contributed_atLeast_ONECE)[1])

write.table(RESULTS_SUMMARY,"~/Documents/PLINK/Tables/summary_IBD_plink.txt",quote=F,row.names=F,col.name=F,sep="\t")




##PLOT
# 1. Divide all segments by landrace population.
# 2. For each chromosomal segment. Identify the wild populations contibuting.Plot them with different colors and at different levels.

#Import the summary table from STRUCTURE analysis (Clumpp_output.R)
Land_str<-read.table("~/Documents/land_assignment_k4_iSelect6135.txt",header=T)

Wild_str<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)

#1. Divide in landrace pops.

Land1<-subset(Land_str, Land_str$Populations == "1")
Land2<-subset(Land_str, Land_str$Populations == "2")
Land3<-subset(Land_str, Land_str$Populations == "3")
Land4<-subset(Land_str, Land_str$Populations == "4")


if (identical(length(grep("^WBDC",All_important_segments[,3])),dim(All_important_segments)[1])) {

all_genetic_pos<-genmap_iSelect[(genmap_iSelect$SNP %in% genmap_30SNPs[,2]),]

#starting plot
GENMAP_plot<-cbind(all_genetic_pos[,6],c(rep(1,dim(all_genetic_pos)[1])),c(rep(2,dim(all_genetic_pos)[1])),c(rep(3,dim(all_genetic_pos)[1])),c(rep(4,dim(all_genetic_pos)[1])))

pdf("~/Documents/PLINK/Plots/IBD_wild_land_plink_30SNPs.pdf",width=11,height=6)
par(mar=c(4,9,4,1),xpd=F)

plot(GENMAP_plot[,1], GENMAP_plot[,c(2)],ylim=c(1,6),xlim=c(1,1113),cex=0.2,col="black",ylab="",xlab="Linkage group",pch=20,xaxt="n",yaxt="n",main="IBS segments between wild and barley landraces, extending 30SNPs ")
points(GENMAP_plot[,1], GENMAP_plot[,3],ylim=c(1,6),cex=0.2,col="black",pch=20)
points(GENMAP_plot[,1], GENMAP_plot[,4],ylim=c(1,6),cex=0.2,col="black",pch=20)
points(GENMAP_plot[,1], GENMAP_plot[,5],ylim=c(1,6),cex=0.2,col="black",pch=20)

Land1_imp_segments<-All_important_segments[(All_important_segments[,2] %in% Land1$Samples),]
Land2_imp_segments<-All_important_segments[(All_important_segments[,2] %in% Land2$Samples),]
Land3_imp_segments<-All_important_segments[(All_important_segments[,2] %in% Land3$Samples),]
Land4_imp_segments<-All_important_segments[(All_important_segments[,2] %in% Land4$Samples),]

wild1<-subset(Wild_str, Wild_str$pop6_hap == '1')
wild2<-subset(Wild_str, Wild_str$pop6_hap == '2')
wild3<-subset(Wild_str, Wild_str$pop6_hap == '3')
wild5<-subset(Wild_str, Wild_str$pop6_hap == '5')
wild6<-subset(Wild_str, Wild_str$pop6_hap == '6')

#For each landrace population at the time
for (l in 1:4){

land<-get(paste("Land",l,"_imp_segments",sep=""))

	#separate by wild population
	land_w1<-land[(land[,3] %in% wild1$Sample),]
	land_w2<-land[(land[,3] %in% wild2$Sample),]
	land_w3<-land[(land[,3] %in% wild3$Sample),]
	land_w5<-land[(land[,3] %in% wild5$Sample),]
	land_w6<-land[(land[,3] %in% wild6$Sample),]
	
	#for each ancestor find the start and end point of every segment contributed
	

		#find unique segments shared: Northen Mesopotamia
		if (dim(land_w1)[1] >0) {
		UNIQUE_SEGMENTS<-unique(sort(land_w1[,1]))
		
		#Make a table of how many wild and how many different landraces interact at each segment for each landrace population
		SHARED_w1<-matrix(ncol=2,nrow=(length(UNIQUE_SEGMENTS)))
		colnames(SHARED_w1) <-c("Landraces","Wild")
		rownames(SHARED_w1)<-as.data.frame(UNIQUE_SEGMENTS)[,1]
		head(SHARED_w1)
		for (s in 1:length(UNIQUE_SEGMENTS)) {
	
		SEG<-UNIQUE_SEGMENTS[s]
		Seg_ind<-grep(SEG,land_w1[,1])
		#which land and wild in this segment?
		land_wild<-land_w1[Seg_ind,]
		SHARED_w1[s,1]<-length(unique(sort(land_wild[,2])))
		SHARED_w1[s,2]<-length(unique(sort(land_wild[,3])))
		}

		SINGLETONS<-unique(sort(c(which(SHARED_w1[,2] <=0) )))
		Frequent_w1<-SHARED_w1[-SINGLETONS,]
		
		UNIQUE_SEGMENTS<-UNIQUE_SEGMENTS[-SINGLETONS]
		if (length(UNIQUE_SEGMENTS) >0) {
		
		#for each segment find the cM position for start and end positions
		for (u in 1:length(UNIQUE_SEGMENTS)) {
			Segment_analyzed<-which(SEGMENT_information[,1] == as.character(UNIQUE_SEGMENTS[u]))
			SNP_start<-SEGMENT_information [Segment_analyzed,4]
			SNP_end<-SEGMENT_information [Segment_analyzed,10]
			genetic_SNP_start<-which(genmap_iSelect$SNP == as.character(SNP_start))
			genetic_SNP_end<-which(genmap_iSelect$SNP == as.character(SNP_end))
			Position_start<-genmap_iSelect[c(genetic_SNP_start),]
			Position_end<-genmap_iSelect[c(genetic_SNP_end),]
			
			segments(Position_start$Cumulative,(l+0.09), Position_end$Cumulative,(l+0.09),col="red")
				}
			}
		}
		

		#find unique segments shared:Southern Levant
		if (dim(land_w2)[1] >0) {
		UNIQUE_SEGMENTS<-unique(sort(land_w2[,1]))
		
		#Make a table of how many wild and how many different landraces interact at each segment for each landrace population
		SHARED_w2<-matrix(ncol=2,nrow=(length(UNIQUE_SEGMENTS)))
		colnames(SHARED_w2) <-c("Landraces","Wild")
		rownames(SHARED_w2)<-as.data.frame(UNIQUE_SEGMENTS)[,1]
		head(SHARED_w2)
		for (s in 1:length(UNIQUE_SEGMENTS)) {
	
		SEG<-UNIQUE_SEGMENTS[s]
		Seg_ind<-grep(SEG,land_w2[,1])
		#which land and wild in this segment?
		land_wild<-land_w2[Seg_ind,]
		SHARED_w2[s,1]<-length(unique(sort(land_wild[,2])))
		SHARED_w2[s,2]<-length(unique(sort(land_wild[,3])))
		}

		SINGLETONS<-unique(sort(c(which(SHARED_w2[,2] <= 0) )))
		Frequent_w2<-SHARED_w2[-SINGLETONS,]
		
		UNIQUE_SEGMENTS<-UNIQUE_SEGMENTS[-SINGLETONS]
		if (length(UNIQUE_SEGMENTS) >0) {
		
		
				
		#for each segment find the cM position for start and end positions
		for (u in 1:length(UNIQUE_SEGMENTS)) {
			Segment_analyzed<-which(SEGMENT_information[,1] == as.character(UNIQUE_SEGMENTS[u]))
			SNP_start<-SEGMENT_information [Segment_analyzed,4]
			SNP_end<-SEGMENT_information [Segment_analyzed,10]
			genetic_SNP_start<-which(genmap_iSelect$SNP == as.character(SNP_start))
			genetic_SNP_end<-which(genmap_iSelect$SNP == as.character(SNP_end))
			Position_start<-genmap_iSelect[c(genetic_SNP_start),]
			Position_end<-genmap_iSelect[c(genetic_SNP_end),]
			segments(Position_start$Cumulative,(l+0.05), Position_end$Cumulative,(l+0.05),col="blue")
				}
			}
		}
		
		#find unique segments shared with Syrian Desert
		if (dim(land_w3)[1] >0) {
		UNIQUE_SEGMENTS<-unique(sort(land_w3[,1]))
		
		#Make a table of how many wild and how many different landraces interact at each segment for each landrace population
		SHARED_w3<-matrix(ncol=2,nrow=(length(UNIQUE_SEGMENTS)))
		colnames(SHARED_w3) <-c("Landraces","Wild")
		rownames(SHARED_w3)<-as.data.frame(UNIQUE_SEGMENTS)[,1]
		head(SHARED_w3)
		for (s in 1:length(UNIQUE_SEGMENTS)) {
	
		SEG<-UNIQUE_SEGMENTS[s]
		Seg_ind<-grep(SEG,land_w3[,1])
		#which land and wild in this segment?
		land_wild<-land_w3[Seg_ind,]
		SHARED_w3[s,1]<-length(unique(sort(land_wild[,2])))
		SHARED_w3[s,2]<-length(unique(sort(land_wild[,3])))
		}

		SINGLETONS<-unique(sort(c(which(SHARED_w3[,2] <=0) )))
		Frequent_w3<-SHARED_w3[-SINGLETONS,]
		
		UNIQUE_SEGMENTS<-UNIQUE_SEGMENTS[-SINGLETONS]
		if (length(UNIQUE_SEGMENTS) >0) {
		
		#for each segment find the cM position for start and end positions
		for (u in 1:length(UNIQUE_SEGMENTS)) {
			Segment_analyzed<-which(SEGMENT_information[,1] == as.character(UNIQUE_SEGMENTS[u]))
			SNP_start<-SEGMENT_information [Segment_analyzed,4]
			SNP_end<-SEGMENT_information [Segment_analyzed,10]
			genetic_SNP_start<-which(genmap_iSelect$SNP == as.character(SNP_start))
			genetic_SNP_end<-which(genmap_iSelect$SNP == as.character(SNP_end))
			Position_start<-genmap_iSelect[c(genetic_SNP_start),]
			Position_end<-genmap_iSelect[c(genetic_SNP_end),]
			segments(Position_start$Cumulative,(l+0.07), Position_end$Cumulative,(l+0.07),col="green")
				}
			}
		}
		
		#find unique segments shared with Northen Levant
		if (dim(land_w5)[1] >0) {
		UNIQUE_SEGMENTS<-unique(sort(land_w5[,1]))
		
		#Make a table of how many wild and how many different landraces interact at each segment for each landrace population
		SHARED_w5<-matrix(ncol=2,nrow=(length(UNIQUE_SEGMENTS)))
		colnames(SHARED_w5) <-c("Landraces","Wild")
		rownames(SHARED_w5)<-as.data.frame(UNIQUE_SEGMENTS)[,1]
		head(SHARED_w5)
		for (s in 1:length(UNIQUE_SEGMENTS)) {
	
		SEG<-UNIQUE_SEGMENTS[s]
		Seg_ind<-grep(SEG,land_w5[,1])
		#which land and wild in this segment?
		land_wild<-land_w5[Seg_ind,]
		SHARED_w5[s,1]<-length(unique(sort(land_wild[,2])))
		SHARED_w5[s,2]<-length(unique(sort(land_wild[,3])))
		}

		SINGLETONS<-unique(sort(c(which(SHARED_w5[,2] <= 0) )))
		Frequent_w5<-SHARED_w5[-SINGLETONS,]
		
		UNIQUE_SEGMENTS<-UNIQUE_SEGMENTS[-SINGLETONS]
		if (length(UNIQUE_SEGMENTS) >0) {
		
		
		#for each segment find the cM position for start and end positions
		for (u in 1:length(UNIQUE_SEGMENTS)) {
			Segment_analyzed<-which(SEGMENT_information[,1] == as.character(UNIQUE_SEGMENTS[u]))
			SNP_start<-SEGMENT_information [Segment_analyzed,4]
			SNP_end<-SEGMENT_information [Segment_analyzed,10]
			genetic_SNP_start<-which(genmap_iSelect$SNP == as.character(SNP_start))
			genetic_SNP_end<-which(genmap_iSelect$SNP == as.character(SNP_end))
			Position_start<-genmap_iSelect[c(genetic_SNP_start),]
			Position_end<-genmap_iSelect[c(genetic_SNP_end),]
			segments(Position_start$Cumulative,(l+0.11), Position_end$Cumulative,(l+0.11),col="orange")
				}
			}
		}
		
		#find unique segments shared with Northen Levant
		if (dim(land_w6)[1] >0) {
		UNIQUE_SEGMENTS<-unique(sort(land_w6[,1]))
		
		#Make a table of how many wild and how many different landraces interact at each segment for each landrace population
		SHARED_w6<-matrix(ncol=2,nrow=(length(UNIQUE_SEGMENTS)))
		colnames(SHARED_w6) <-c("Landraces","Wild")
		rownames(SHARED_w6)<-as.data.frame(UNIQUE_SEGMENTS)[,1]
		head(SHARED_w6)
		for (s in 1:length(UNIQUE_SEGMENTS)) {
	
		SEG<-UNIQUE_SEGMENTS[s]
		Seg_ind<-grep(SEG,land_w6[,1])
		#which land and wild in this segment?
		land_wild<-land_w6[Seg_ind,]
		SHARED_w6[s,1]<-length(unique(sort(land_wild[,2])))
		SHARED_w6[s,2]<-length(unique(sort(land_wild[,3])))
		}

		SINGLETONS<-unique(sort(c(which(SHARED_w6[,2] <= 0) )))
		Frequent_w6<-SHARED_w6[-SINGLETONS,]
		UNIQUE_SEGMENTS<-UNIQUE_SEGMENTS[-SINGLETONS]
		if (length(UNIQUE_SEGMENTS) >0) {

		
		#for each segment find the cM position for start and end positions
		for (u in 1:length(UNIQUE_SEGMENTS)) {
			Segment_analyzed<-which(SEGMENT_information[,1] == as.character(UNIQUE_SEGMENTS[u]))
			SNP_start<-SEGMENT_information [Segment_analyzed,4]
			SNP_end<-SEGMENT_information [Segment_analyzed,10]
			genetic_SNP_start<-which(genmap_iSelect$SNP == as.character(SNP_start))
			genetic_SNP_end<-which(genmap_iSelect$SNP == as.character(SNP_end))
			Position_start<-genmap_iSelect[c(genetic_SNP_start),]
			Position_end<-genmap_iSelect[c(genetic_SNP_end),]
			segments(Position_start$Cumulative,(l+0.03), Position_end$Cumulative,(l+0.03),col="purple")
				}
			}
		}
		
	IBD_LONG_freq<-rbind(Frequent_w1, Frequent_w2, Frequent_w3, Frequent_w5, Frequent_w6)
	assign(paste("Land",l,"_IBD_long_freq",sep=""), IBD_LONG_freq)
}

abline(v=c(146.9, 325.23, 489.65, 620.55, 804.38,943.77),lty=4,col="gray")

legend("topleft", c("Northern Levant","Northern Mesopotamia","Syrian Desert","Southern Levant","Central Asia", "Linkage group boundaries","SNP"),lwd=2, cex=0.70 , col=(c("orange","red","green","blue","purple","gray","black")),lty=c(rep(1,5),5,NA),pch=c(NA,NA,NA,NA,NA,NA,20),bg="white",title=expression(bold("Wild Populations")))

axis(side=1,at=c(70,230,400,555,710,870,1050),c("1H","2H","3H","4H","5H","6H","7H"))
axis(side=2,at=c(1,2,3,4),c("Central European","Asia","Coastal Mediterranean","East Africa"),cex.axis=0.9,las=1)
#dev.off()

} else print("ERROR:Landraces and Wild samples are mixed in one column")

