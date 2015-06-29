
#Make tables for each chromosome. This has the solved heterozygotes calls when the prob of making the right call was >90%. All other sites were considered as NA, and what was missing data was replaced for NA. 

#By hand separate the results in file My_output_1 (where 1 is the chromosome number):
#	⁃	separate in files 1) the inferred genotypes BESTPAIRS1->(Chr1_phasedHaps.txt) by hand replace the 0 # for just # 
#	⁃	2)the table with === ? and () [] , called BESTPAIRS2. (Chr1_originalPhase.txt)Also replace the 0 # for just #
	
rm(list=ls())
#Select here the chromosome to work with:
for (CHR in 1:7){

DATA <-read.table(paste("~/PHASE/output/Separate_output/Chr",CHR,"_originalPhase.txt",sep=""),header=F)

dim(DATA)


PHASED<-read.table(paste("~/PHASE/output/Chr",CHR,"_phasedHaps.txt",sep=""),header=F)


dim(PHASED)
PHASED->PHASED_DATA

#Bring SNP names per chromosome
INPUT_phase_chr<-read.table(paste("~/Documents/land_prephase_chr_",CHR,".txt",sep=""),header=T,row.names=1)
#INPUT_phase_chr<-read.table(paste("~/Documents/wild_prephase_chr_",CHR,".txt",sep=""),header=T,row.names=1)

SNP_chr<-colnames(INPUT_phase_chr)
Samples<-row.names(INPUT_phase_chr)

for (c in 1:(dim(PHASED_DATA)[2])) {
	for (s in 1:(dim(PHASED_DATA)[1])) {
		#probabilities less than 90% for infered (difficult to resolve hete) values are marked with "( )", we will replace these values with Missing data = NA
	if (DATA[s,c] == '?') (PHASED_DATA[s,c] <- NA) else if (as.data.frame(DATA[s,c]) == "(A)") (PHASED_DATA[s,c] <- NA) else if (as.data.frame(DATA[s,c]) == "(B)") (PHASED_DATA[s,c] <- NA)
	}
}

colnames(PHASED_DATA)<-SNP_chr
row.names(PHASED_DATA)<-Samples
PHASED_DATA[300:310,1:20]

write.table(PHASED_DATA,paste("~/PHASE/output/Separate_output/PHASED_hete_land/Chr",CHR,"_land_noHete.txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
#write.table(PHASED_DATA,paste("~/PHASE/output/Separate_output/PHASED_hete_wild/Chr",CHR,"_wild_noHete.txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")



#==============Put all chromosomes together==============

if (FALSE){
	
	CHR1<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Chr1_land_noHete.txt",header=T,row.names=1)
	CHR2<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Chr2_land_noHete.txt",header=T,row.names=1)
	CHR3<-read.table("~/PHASE/outputSeparate_output/PHASED_hete_land/Chr3_land_noHete.txt",header=T,row.names=1)
	CHR4<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Chr4_land_noHete.txt",header=T,row.names=1)
	CHR5<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Chr5_land_noHete.txt",header=T,row.names=1)
	CHR6<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Chr6_land_noHete.txt",header=T,row.names=1)
	CHR7<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Chr7_land_noHete.txt",header=T,row.names=1)
	
	ALL_chr<-	cbind(CHR1,CHR2,CHR3,CHR4,CHR5,CHR6,CHR7)
	write.table(ALL_chr,"~/PHASE/output/Separate_output/PHASED_hete_land/Phase_diploid_land.txt",quote=F,row.names=T,col.names=T,sep="\t")
}

if (FALSE){
	
	CHR1<-read.table("~/PHASE/output/Separate_output/PHASED_hete_wild/Chr1_wild_noHete.txt",header=T,row.names=1)
	CHR2<-read.table("~/PHASE/output/output/Separate_output/PHASED_hete_wild/Chr2_wild_noHete.txt",header=T,row.names=1)
	CHR3<-read.table("~/PHASE/output/output/Separate_output/PHASED_hete_wild/Chr3_wild_noHete.txt",header=T,row.names=1)
	CHR4<-read.table("~/PHASE/output/Separate_output/PHASED_hete_wild/Chr4_wild_noHete.txt",header=T,row.names=1)
	CHR5<-read.table("~/PHASE/output/Separate_output/PHASED_hete_wild/Chr5_wild_noHete.txt",header=T,row.names=1)
	CHR6<-read.table("~/PHASE/output/output/Separate_output/PHASED_hete_wild/Chr6_wild_noHete.txt",header=T,row.names=1)
	CHR7<-read.table("~/PHASE/output/output/Separate_output/PHASED_hete_wild/Chr7_wild_noHete.txt",header=T,row.names=1)
	
	ALL_chr<-	cbind(CHR1,CHR2,CHR3,CHR4,CHR5,CHR6,CHR7)
	write.table(ALL_chr,"~/PHASE/output/Separate_output/PHASED_hete_wild/Phase_diploid_wild.txt",quote=F,row.names=T,col.names=T,sep="\t")
}
