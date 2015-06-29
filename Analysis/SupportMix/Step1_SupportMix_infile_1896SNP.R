#This scripts will transform PHASED, no IMPUTED, files to input for SupportMix. For 1896SNPs common between wild and landraces.
#Run the code. the files with be sent to pre-created hapmap2Subset and tpedData folders here:~/Dropbox/ANITA_301213/THESIS/Data/for_SupporMix/Supportmix_barley_phased_1896/

rm(list=ls())

#1. Import the data for each chromosome and put all chrs together for all landraces, then for each wild population.
#2. Turn the table to have SNPs in rows and samples in columns
#3. Identify the SNPs in the iSelect map to get ther positions.
#4. Create formated tped files and write them to the right folder
#5. Create tfam files
#6. Create tped and tfam files for each wild populaton separately 

Chr1_land <-read.table("~/PHASE/PHASED_hete_land/Chr1_land_noHete.txt",header=T,row.names=1)
Chr2_land <-read.table("~/PHASE/PHASED_hete_land/Chr2_land_noHete.txt",header=T,row.names=1)
Chr3_land <-read.table("~/PHASE/PHASED_hete_land/Chr3_land_noHete.txt",header=T,row.names=1)
Chr4_land <-read.table("~/PHASE/PHASED_hete_land/Chr4_land_noHete.txt",header=T,row.names=1)
Chr5_land <-read.table("~/PHASE/PHASED_hete_land/Chr5_land_noHete.txt",header=T,row.names=1)
Chr6_land <-read.table("~/PHASE/PHASED_hete_land/Chr6_land_noHete.txt",header=T,row.names=1)
Chr7_land <-read.table("~/PHASE/PHASED_hete_land/Chr7_land_noHete.txt",header=T,row.names=1)

LAND_phased<-cbind(Chr1_land, Chr2_land, Chr3_land, Chr4_land, Chr5_land, Chr6_land, Chr7_land)
LAND_phased[1:10,1:20]

#Replace NA with 0 =  plink format for missing data, doing it by rows = 1, it will turn the table to have samples in cols and SNPs in rows
replaceNA<-function(dat){
	dat[is.na(dat)]<-'0'
	return(dat)
	}
LAND_phased_t<-as.data.frame(apply(LAND_phased,1,replaceNA))
dim(LAND_phased_t)

#Call genetic map and select SNPs that are in the LANDraces

genmap<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)
head(genmap)

genmap_land<-genmap[(genmap$SNP %in% row.names(LAND_phased_t)),]
dim(genmap_land)

#Format tped file
extra<-rep(0,(dim(LAND_phased_t)[1]))
land_admix_tped<-cbind(as.data.frame(genmap_land$Chromosome),as.data.frame(genmap_land$SNP),extra,as.data.frame(genmap_land$cm),as.data.frame(LAND_phased_t))

write.table(land_admix_tped,"~/Documents/Supportmix_barley_phased_1896/input_HEADERS/land_admixed_k5.tped",quote=F,row.names=F,col.names=T,sep=" ")
write.table(land_admix_tped,"~/Documents/Supportmix_barley_phased_1896/tpedData/land_admixed_k5.tped",quote=F,row.names=F,col.names=F,sep=" ")

####tfam file for landraces

SAMPLES_land<-as.data.frame(colnames(LAND_phased_t))
#Select odd numbers, so we have only one name for each sample.

Odd_sample<-seq(1,nrow(LAND_phased), by=2)

Sample_name_land<-SAMPLES_land[Odd_sample,]

extra_col<-rep(0,length(Sample_name_land))
land_admix_tfam<-cbind(as.data.frame(Sample_name_land), as.data.frame(Sample_name_land),extra_col, extra_col, extra_col, extra_col)

write.table(land_admix_tfam,"~/Documents/Supportmix_barley_phased_1896/tpedData/land_admixed_k5.tfam",quote=F,row.names=F,col.names=F,sep="\t")

######make a genetic map file for each chromosome
genetic_files<-cbind(genmap_land[,4],rep(1,(dim(genmap_land)[1])), genmap_land[,4], genmap_land[,3])
#Separate by chromosome and remove last column
chr1_genmap<-subset(genetic_files, genetic_files[,4] == '1')[,c(1:3)]
chr2_genmap<-subset(genetic_files, genetic_files[,4] == '2')[,c(1:3)]
chr3_genmap<-subset(genetic_files, genetic_files[,4] == '3')[,c(1:3)]
chr4_genmap<-subset(genetic_files, genetic_files[,4] == '4')[,c(1:3)]
chr5_genmap<-subset(genetic_files, genetic_files[,4] == '5')[,c(1:3)]
chr6_genmap<-subset(genetic_files, genetic_files[,4] == '6')[,c(1:3)]
chr7_genmap<-subset(genetic_files, genetic_files[,4] == '7')[,c(1:3)]

for (i in 1:7){
INPUT<-get(paste("chr",i,"_genmap",sep=""))
OUTPUT<-paste("~/Documents/Supportmix_barley_phased_1896/hapmap2Subset/genetic_map_chr",i,"_file.txt",sep="")
write.table(INPUT,OUTPUT,quote=F,row.names=F,col.names=F,sep="\t")
	}
##################----------WILD TPED Tfam FILES-------------------------################
#1. Get phased data for all chr
#2. Replace NA for '0' and turn the table
#3. Select markers in genmap and make tmax tables 


Chr1_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr1_wild_noHete.txt",header=T,row.names=1)
Chr2_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr2_wild_noHete.txt",header=T,row.names=1)
Chr3_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr3_wild_noHete.txt",header=T,row.names=1)
Chr4_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr4_wild_noHete.txt",header=T,row.names=1)
Chr5_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr5_wild_noHete.txt",header=T,row.names=1)
Chr6_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr6_wild_noHete.txt",header=T,row.names=1)
Chr7_wild <-read.table("~/PHASE/PHASED_hete_wild/Chr7_wild_noHete.txt",header=T,row.names=1)

WILD_phased<-cbind(Chr1_wild, Chr2_wild, Chr3_wild, Chr4_wild, Chr5_wild, Chr6_wild, Chr7_wild)
WILD_phased[1:10,1:20]
dim(WILD_phased)

#Replace NA for 0 and turn table to have SNPs in rows and samples names in columns
WILD_phased_t<-as.data.frame(apply(WILD_phased,1,replaceNA))
dim(WILD_phased_t)

#Call genetic map and select SNPs that are in the LANDraces

genmap_WILD<-genmap[(genmap$SNP %in% row.names(WILD_phased_t)),]
dim(genmap_WILD)

#Get columns for SNP information
extra<-rep(0,(dim(WILD_phased_t)[1]))
WILD_SNP_info<-cbind(as.data.frame(genmap_WILD$Chromosome),as.data.frame(genmap_WILD$SNP),extra,as.data.frame(genmap_WILD$cm))

#Separate wild barley by populations
str_wild<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)
head(str_wild)

#make another str_wild file adding the '_2' to each sample name. Then combine the files to later separate all the wild samples.
str_wild_2<-cbind(paste(str_wild[,1],'_2',sep=""),as.data.frame(str_wild[,2:4]))
colnames(str_wild_2)<-c("Sample","pop6_hap","latitude","longitude")

wild_str_all<-rbind(as.data.frame(str_wild),as.data.frame(str_wild_2))

#Separate wild sample populations and add SNP information
pop1<-subset(wild_str_all, wild_str_all$pop6_hap =='1')
pop2<-subset(wild_str_all, wild_str_all$pop6_hap =='2')
pop3<-subset(wild_str_all, wild_str_all$pop6_hap =='3')
pop5<-subset(wild_str_all, wild_str_all$pop6_hap =='5')
pop6<-subset(wild_str_all, wild_str_all$pop6_hap =='6')

Wild1_phased<-cbind(WILD_SNP_info ,as.data.frame(WILD_phased_t[,(colnames(WILD_phased_t) %in% pop1[,1])]))
Wild2_phased<-cbind(WILD_SNP_info ,as.data.frame(WILD_phased_t[,(colnames(WILD_phased_t) %in% pop2[,1])]))
Wild3_phased<-cbind(WILD_SNP_info ,as.data.frame(WILD_phased_t[,(colnames(WILD_phased_t) %in% pop3[,1])]))
Wild5_phased<-cbind(WILD_SNP_info ,as.data.frame(WILD_phased_t[,(colnames(WILD_phased_t) %in% pop5[,1])]))
Wild6_phased<-cbind(WILD_SNP_info ,as.data.frame(WILD_phased_t[,(colnames(WILD_phased_t) %in% pop6[,1])]))


#Although we generate the input files for the Caspian Sea population, this is not taking into account when running SupportMix
for (i in c(1,2,3,5,6)){
	INPUT<-get(paste("Wild",i,"_phased",sep=""))
	OUTPUT<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/tpedData/ancestral_pop",i,".tped",sep="")
	Output_header<-paste("~/Documentes/SupporMix/Supportmix_barley_phased_1896/input_HEADERS/ancestral_pop",i,".tped",sep="")
	write.table(INPUT,OUTPUT,quote=F,row.names=F,col.names=F,sep=" ")
	write.table(INPUT, Output_header,quote=F,row.names=F,col.names=T,sep=" ")
}

#### tfam for WILD
#for each wild population, take the Wildx_phase file and grab the samples names. Create an extra column of Zeros.
for (i in c(1,2,3,5,6)) {
phased_WILD<-get(paste("Wild",i,"_phased",sep=""))
SAMPLES_wild<-as.data.frame(as.data.frame(colnames(phased_WILD))[-c(1:4),1])
#Select odd numbers, so we have only one name for each sample.
Odd_sample_wild<-seq(1,nrow(SAMPLES_wild), by=2)

Sample_name_wild<-SAMPLES_wild[Odd_sample_wild,1]

extra_col_wild<-rep(0,length(Sample_name_wild))
wild_admix_tfam<-cbind(as.data.frame(Sample_name_wild), as.data.frame(Sample_name_wild), extra_col_wild, extra_col_wild, extra_col_wild, extra_col_wild)

OUTPUT<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/tpedData/ancestral_pop",i,".tfam",sep="")
write.table(wild_admix_tfam,OUTPUT,quote=F,row.name=F,col.name=F,sep="\t")
}


#####go to my home directory where all the tped and ped files are and gunzip them.
# $tpedData gonzales$ gzip *
#Use the config file to run SupportMix
#$ SupportMix -C Config_1896_wildk5.cfg
#Then continue with Step2_Suppormix_output_1896.R