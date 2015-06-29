#This script prepared the input files for PLINK to estimate Identity by Descent using 30 SNPs haplotypes. 
rm(list=ls())

land<-read.table("~/PHASE/output/Separate_output/PHASED_hete_land/Phase_diploid_land.txt",header=T,row.names=1)
wild_all<-read.table("~/PHASE/output/Separate_output/PHASED_hete_wild/Phase_diploid_wild.txt",header=T,row.names=1)

#Remove Caspian Sea wild samples
wild_pops<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)
caspian<-subset(wild_pops, wild_pops$pop6_hap == '4')
wild<-wild_all[!(row.names(wild_all) %in% caspian$Sample),]


barley2_ordered<-rbind(as.data.frame(land),as.data.frame(wild))


replaceNA<-function(dat){
	dat[is.na(dat)]<-'0'
	return(dat)
	}
genotypes<-as.data.frame(apply(barley2_ordered,2,replaceNA))

genotypes[1:10,1:10]



#Separate the two haplotypes

hap1<- seq(1,nrow(genotypes),by=2)
haplotype1<- genotypes[hap1,]

hap2<- seq(2,nrow(genotypes),by=2)
haplotype2<- genotypes[hap2,]

#combine haplotypes to apper SNP1_allele1 SNP1_allele2

removed<-cbind(haplotype1,haplotype2)
removed[1:8,1:10]
dim(removed)
removed_or<-removed[,order(colnames(removed))]
dim(removed_or)
removed_or[1:8,1:10]

#Sort SNPs by genetic position:Turn the table, sort the SNPs based on genetic map, and turn the table again to have SNPs in columns and Samples in rows
t_removed_or<-as.data.frame(t(removed_or))

genmap<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)
head(genmap)

#select 1896 SNPs

genmap_part<-genmap[(genmap$SNP %in% row.names(t_removed_or)),]

Index_hap_odd<-seq(1,(1896*2),2)
genmap_odd<-cbind(genmap_part,Index_hap_odd)
colnames(genmap_odd)<-c("SNP"  , "SNP_2", "Chro", "cm", "index", "Cumulative", "Chromosome" ,"Index_hap")

Index_hap_even<-seq(2,(1896*2),2)
genmap_even<-cbind(paste(genmap_part[,1],".1",sep=""),as.data.frame(genmap_part[,2:7]),Index_hap_even)
colnames(genmap_even)<-c("SNP"  , "SNP_2", "Chro", "cm", "index", "Cumulative", "Chromosome" ,"Index_hap")

genmap_all<-rbind(genmap_odd,genmap_even)
genmap_all_in<-genmap_all[order(genmap_all$Index_hap),]

rownames(genmap_all_in)<-genmap_all_in$SNP

genmap_all_in_sort<-genmap_all_in[order(row.names(genmap_all_in)),]
identical(rownames(t_removed_or),rownames(genmap_all_in_sort))
genotypes_SNP<-cbind(as.data.frame(genmap_all_in_sort), as.data.frame(t_removed_or))
rownames(genotypes_SNP)<-row.names(t_removed_or)

#sort SNPs by Index_haplotype, this column has the right genetic order for both haplotypes
genotypes_SNPor<-genotypes_SNP[order(genotypes_SNP$Index_hap),]
#Remove extra columns for SNP ordering and turn the table to have SNPs in columns

genotypes_SNPor<-genotypes_SNPor[,-c(1:8)]

genotypes_SNPor_t<-as.data.frame(t(genotypes_SNPor))

#Create extra columns for PED file
ZERO<-rep("0",(dim(genotypes_SNPor_t)[1]))
NINE<-rep("-9",(dim(genotypes_SNPor_t)[1]))

PED<-cbind(as.data.frame(ZERO),as.data.frame(row.names(genotypes_SNPor_t)),as.data.frame(ZERO),as.data.frame(ZERO),as.data.frame(ZERO),as.data.frame(NINE),as.data.frame(genotypes_SNPor_t))
PED[1:8,1:10]

#write.table(PED,"~/Documents/Anita/THESIS/Data/for_IBS/Barley_1896/Input/WILD_LAND_header.ped",quote=F,row.names=F,col.names=T,sep="\t")


#Separate by chromosomes
Ped_initial<-PED[,c(1:6)]

snp_single<-cbind(as.data.frame(genmap_part[,c(7,1,4)]),round(genmap_part[,4],digits=0))
head(snp_single)
colnames(snp_single)<-c("Chromosome","SNP","cm","rounded")
snp_double<-cbind(as.data.frame(snp_single[,1]),paste(snp_single[,2],".1",sep=""),as.data.frame(snp_single[,3:4]))
colnames(snp_double)<-c("Chromosome","SNP","cm","rounded")

snp<-rbind(as.data.frame(snp_single),as.data.frame(snp_double))
chr1<-subset(snp,snp[,1] == '1')
chr2<-subset(snp,snp[,1] == '2')
chr3<-subset(snp,snp[,1] == '3')
chr4<-subset(snp,snp[,1] == '4')
chr5<-subset(snp,snp[,1] == '5')
chr6<-subset(snp,snp[,1] == '6')
chr7<-subset(snp,snp[,1] == '7')

PED_chr1<-PED[,(colnames(PED) %in% chr1[,2])]
PED_chr2<-PED[,(colnames(PED) %in% chr2[,2])]
PED_chr3<-PED[,(colnames(PED) %in% chr3[,2])]
PED_chr4<-PED[,(colnames(PED) %in% chr4[,2])]
PED_chr5<-PED[,(colnames(PED) %in% chr5[,2])]
PED_chr6<-PED[,(colnames(PED) %in% chr6[,2])]
PED_chr7<-PED[,(colnames(PED) %in% chr7[,2])]

#Add the initial columns
PED_chr1_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr1))
PED_chr2_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr2))
PED_chr3_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr3))
PED_chr4_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr4))
PED_chr5_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr5))
PED_chr6_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr6))
PED_chr7_ready<-cbind(as.data.frame(Ped_initial),as.data.frame(PED_chr7))

for (i in 1:7){
INPUT<-get(paste("PED_chr",i,"_ready",sep=""))
OUTPUT<-paste("~/Documents/Input/ByChromosome/Wild_land_chr",i,".ped",sep="")
write.table(INPUT,OUTPUT,quote=F,row.names=F,col.names=F,sep="\t")
}

#====For IBD using plink only====

#Devide .map in chromosomes
mapfile<-snp_single
head(mapfile)

map_chr1<-subset(mapfile,mapfile[,1]=='1')
map_chr2<-subset(mapfile,mapfile[,1]=='2')
map_chr3<-subset(mapfile,mapfile[,1]=='3')
map_chr4<-subset(mapfile,mapfile[,1]=='4')
map_chr5<-subset(mapfile,mapfile[,1]=='5')
map_chr6<-subset(mapfile,mapfile[,1]=='6')
map_chr7<-subset(mapfile,mapfile[,1]=='7')

for (i in 1:7){
INPUT<-get(paste("map_chr",i,sep=""))
OUTPUT<-paste("~/Documents/Input/ByChromosome/chr",i,".map",sep="")
write.table(INPUT,OUTPUT,quote=F,row.names=F,col.names=F,sep="\t")
}

##fOR THE LANDRACES by Populations
#Separate the landraces in the for populations, using the output from the STRUCTURE analysis.
str<-read.table("~/Documents/land_assignment_k4_iSelect6135.txt",header=T)
head(str)

pop1<-subset(str,str$Populations == 1)
pop2<-subset(str,str$Populations == 2)
pop3<-subset(str,str$Populations == 3)
pop4<-subset(str,str$Populations == 4)

PED1<-PED[(row.names(PED) %in% pop1$Samples),]
PED2<-PED[(row.names(PED) %in% pop2$Samples),]
PED3<-PED[(row.names(PED) %in% pop3$Samples),]
PED4<-PED[(row.names(PED) %in% pop4$Samples),]

for (p in 1:4){
INPUT<-get(paste("PED",p,sep=""))
OUTPUT<-paste("~/Documents/Input/GenomeWide/Pop",p,"_land.ped",sep="")
write.table(INPUT,OUTPUT,quote=F,row.names=F,col.names=F,sep="\t")
}

#Now use Step2_Plink_IBD_1896_input_map_ped30SNP.R to separate the files in SNP windows for plink