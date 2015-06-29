#This script will calculate the allele frequency of each private allele in each of the landraces. The identification of allelic contribution in the landraces was done using SharedPoly_land_carring_wild_PrivateAlleles_1896.R

rm(list=ls())

Ancestry<-read.table("~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/Private_PresentInLand.txt",header=T,row.names=1)
Ancestry[1:10,1:10]
#Order Private Snps according to genetic map iSelect
Private<-read.table("~/Documents/github/BarleyLandraces/Datasets/Private_allele_ALLpops.txt",header=T)
genmap<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)
genmap_or<-genmap[order(genmap$SNP),]

#Select SNPs in genmap that were private in the wild. Identify which SNPs correspond to which wild Pop

genmap_or_priv<-genmap_or[(genmap_or$SNP %in% Private$SNP),]
Private_or<-Private[order(Private$SNP),]

#Check that the same SNPs in same order are being used

if (identical(as.character(Private_or[,"SNP"]), as.character(genmap_or_priv[,"SNP"]))){

GeneticOrder_PrivateSNPs<-cbind(as.data.frame(Private_or),as.data.frame(genmap_or_priv))
GeneticOrder_PrivateSNPs<-GeneticOrder_PrivateSNPs[,-c(4,5)]
} else print ("Error! SNPs names or order differ")

#1. Devide the landraces into the four populations.
#2. Calculate the frequency of each SNP in each land population
#3. Order the SNPs according to a genetic map iSelect
#4. Plot frequency vs genetic map. Color code by wild population. using code Plot_shradpoly_alleleFreq_1896

#1. Devide the landraces into the four populations

str<-read.table("~/Documents/land_assignment_k4_iSelect6135.txt",header=T,row.names=1)
str[1:10,]

Pop1_land<-subset(str,str$Populations == 1)
Pop2_land<-subset(str,str$Populations == 2)
Pop3_land<-subset(str,str$Populations == 3)
Pop4_land<-subset(str,str$Populations == 4)

Pop1_ancestry<-Ancestry[(row.names(Ancestry) %in% row.names(Pop1_land)),]
Pop2_ancestry<-Ancestry[(row.names(Ancestry) %in% row.names(Pop2_land)),]
Pop3_ancestry<-Ancestry[(row.names(Ancestry) %in% row.names(Pop3_land)),]
Pop4_ancestry<-Ancestry[(row.names(Ancestry) %in% row.names(Pop4_land)),]


#2. Calculate the frequency of each SNP in each land population. Order by SNP name.


POPULATION_LAND<-c("Pop1_ancestry","Pop2_ancestry","Pop3_ancestry","Pop4_ancestry")

for (p in 1:length(POPULATION_LAND)){
POPULATION<-get(POPULATION_LAND[p])
Total_ind<-(dim(POPULATION)[1])*2

frequency_private<-NULL
#For each SNP private = 181
for (i in 1:(dim(POPULATION)[2])) {
	
		#See how many are counts comes out of hetes. These will reduce the count. E.g. there are 50 diploid individuals= 100 hap total. If 10 individuals have Minor allele count in homozygote state= 10*2 (total MAF from diploid homo). If two are hete, count= (count_hap*2)+hete == (10*2)+2. Then MAF = count/total_ind=((10*2)+2)/(50*2)
	
	
	Table_class<-matrix(0,ncol=3,nrow=(dim(POPULATION)[2]))
	colnames(Table_class)<-c("major","minor","hete")
	row.names(Table_class)<-colnames(POPULATION)

	#for each SNP i count MAF and order if it comes from homozygote or hererozygote individual. SNP are in col, samples in Rows.
	TABLE_COUNT<-as.data.frame(table(POPULATION[,i]))
	#for item in table order it
	for (t in 1:(dim(TABLE_COUNT)[1])){
		Table_class[i,t] <-TABLE_COUNT[t,2]
	}
	
	#Sum the total number of alleles from the wild./total sample size (=2N) 
	frequency_private[i]<-(Table_class[i,2]*2 + Table_class[i,3])/Total_ind
	
}


frequency_private<-as.data.frame(frequency_private)

frequency_private$SNP_name<-colnames(POPULATION)
frequency_private<-frequency_private[,c(2,1)]
frequency_private_or<-frequency_private[order(frequency_private[,1]),]

#3. Order the SNPs according to a genetic map iSelect in GeneticOrder_PrivateSNPs
if (identical(as.character(frequency_private_or[,"SNP_name"]), as.character(GeneticOrder_PrivateSNPs[,"SNP"]))){
Frequency_genmap<-cbind(as.data.frame(frequency_private_or),as.data.frame(GeneticOrder_PrivateSNPs))

Frequency_genmap_or<-Frequency_genmap[order(Frequency_genmap$index),]

assign(paste("LAND_",p,"_private_freq",sep=""), Frequency_genmap_or)
INPUT<-get(paste("LAND_",p,"_private_freq",sep=""))
OUTPUT<-paste("~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/","LAND_",p,"_private_freq.txt",sep="")
write.table(INPUT,OUTPUT,quote=F,row.names=F,col.names=T,sep="\t")
	} else print ("ERROR! Snps in frequency table are in different order than the genetic map ordered by SNP name table")
}

#Make table freq_private_inEachLand.xls

All_LAND_private_freq_all_info <-cbind(LAND_1_private_freq, LAND_2_private_freq, LAND_3_private_freq, LAND_4_private_freq)
All_land_private_freq<-cbind(as.data.frame(All_LAND_private_freq_all_info[,c(6,1,5,2,12,22,32)]))

colnames(All_land_private_freq)<-c("LinkageGroup","PrivateAllele","WILDpop","CoastalEuropean","Asia","CoastalMediterranean","EastAfrica")

#Separate by wild pop private alleles.
Wild_pop_private<-table(All_land_private_freq$WILDpop)
N_LEVANT<-subset(All_land_private_freq ,All_land_private_freq$WILDpop == "N_Levant")
N_MESOPOTAMIA<-subset(All_land_private_freq ,All_land_private_freq$WILDpop == "N_Mesopotamia")
S_DESERT<-subset(All_land_private_freq ,All_land_private_freq$WILDpop == "S_Desert")
S_LEVANT<-subset(All_land_private_freq ,All_land_private_freq$WILDpop == "S_Levant")
C_ASIA<-subset(All_land_private_freq ,All_land_private_freq$WILDpop == "C_Asia")

#To separate pops in the table
n_lev_sep<-t(as.data.frame(c("Northern Levant","","","","","","")));colnames(n_lev_sep)<-c("LinkageGroup","PrivateAllele","WILDpop","CoastalEuropean","Asia","CoastalMediterranean","EastAfrica")
n_meso_sep<-t(as.data.frame(c("Northern Mesopotamia","","","","","","")));colnames(n_meso_sep)<-c("LinkageGroup","PrivateAllele","WILDpop","CoastalEuropean","Asia","CoastalMediterranean","EastAfrica")
s_desert_sep<-t(as.data.frame(c("Syrian Desert","","","","","","")));colnames(s_desert_sep)<-c("LinkageGroup","PrivateAllele","WILDpop","CoastalEuropean","Asia","CoastalMediterranean","EastAfrica")
s_lev_sep<-t(as.data.frame(c("Southern Levant","","","","","","")));colnames(s_lev_sep)<-c("LinkageGroup","PrivateAllele","WILDpop","CoastalEuropean","Asia","CoastalMediterranean","EastAfrica")
c_asia_sep<-t(as.data.frame(c("Central Asia","","","","","","")));colnames(c_asia_sep)<-c("LinkageGroup","PrivateAllele","WILDpop","CoastalEuropean","Asia","CoastalMediterranean","EastAfrica")

#Put the table together

TABLE_FREQ<-rbind(n_lev_sep,as.data.frame(as.matrix(N_LEVANT)),n_meso_sep,as.data.frame(as.matrix(N_MESOPOTAMIA)),s_desert_sep,as.data.frame(as.matrix(S_DESERT)),s_lev_sep,as.data.frame(as.matrix(S_LEVANT)),c_asia_sep,as.data.frame(as.matrix(C_ASIA)))

#Remove WILDpop column
TABLE_FREQ<-TABLE_FREQ[,-c(3)]

write.table(TABLE_FREQ,"~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/Table_freq_private_inEachLand.xls",quote=F,row.names=F,col.names=T,sep="\t")
#Continue for plot with Plot_shradpoly_alleleFreq_1896.R