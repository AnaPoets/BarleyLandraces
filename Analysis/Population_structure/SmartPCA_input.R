#Create input files for smpartPCA. Using the ancestrymap format:
# Use minor allele as the reference allele
# Set gender to "U" Unknown 
# Set population group as 1 for all, to have an independent result from Structure. We don't have case/controls
# The genotype file uses "9" for missing data and "0" if reference allele (minor allele) is not present


rm(list=ls())
genotype <-read.table("~/Documents/github/BarleyLandraces/Datasets/Land_6152_SNPs_AB.txt",header=T,row.names=1)

dim(genotype)

#Sort SNPs by genetic position

genmap<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)
head(genmap)

genmap_or<-genmap[order(genmap$SNP),]

#Select SNPs that are in genmap and order them

genotypes_map<-genotype[,(colnames(genotype) %in% genmap_or$SNP)]
dim(genotypes_map)

genotypes_t<-as.data.frame(t(genotypes_map))
genotypes_or<-genotypes_t[order(row.names(genotypes_t)),]

#find genmap SNPs that are in the genotypes.And orther by cM position using the index.
genmap_part<-genmap_or[(genmap_or$SNP %in% row.names(genotypes_or)),]

if (identical(as.character(genmap_part$SNP), as.character(row.names(genotypes_or))) == FALSE ) stop (print ("ERROR! SNPs are in different order"))
genotypes_genmap<-cbind(as.data.frame(genmap_part),as.data.frame(genotypes_or))
genotypes_genmap_or<-genotypes_genmap[order(genotypes_genmap$index),]

#Find Minor allele at each SNP
minorAllele<-function(dat){
	AlleleA<-length(which(dat == 'AA'))
	AlleleB<-length(which(dat == 'BB'))
	minor_allele<-if(AlleleA <=AlleleB) {'A'} else {'B'}
	return(minor_allele)
}

MinorAllele_ref<-apply(genotypes_genmap_or,1, minorAllele)
MinorAllele_ref<-as.data.frame(MinorAllele_ref)

#Make a list of Major allele
MajorAllele<-function(dat){
	major_allele <-if (dat == 'A') {'B'} else {'A'}
	return (major_allele)
	}
Major_Allele_ref<-apply(MinorAllele_ref,1, MajorAllele)


dim(genotypes_genmap_or)

#Remove extra columns SNP_2
GENOTYPE_READY<-cbind(as.data.frame(genotypes_genmap_or[,c(1,3,4,5,6)]),as.data.frame(MinorAllele_ref),as.data.frame(genotypes_genmap_or[,8:dim(genotypes_genmap_or)[2]]))

##1. Create SNP.snp file
#Since we don't have physical map, I will create a list of consecutive number for each chromosome.
chr1<-length(which(GENOTYPE_READY$Chro == '1H'))
chr2<-length(which(GENOTYPE_READY$Chro == '2H'))
chr3<-length(which(GENOTYPE_READY$Chro == '3H'))
chr4<-length(which(GENOTYPE_READY$Chro == '4H'))
chr5<-length(which(GENOTYPE_READY$Chro == '5H'))
chr6<-length(which(GENOTYPE_READY$Chro == '6H'))
chr7<-length(which(GENOTYPE_READY$Chro == '7H'))

#make sure that All the chromosomes are present at last once in  GENOTYPE_READY. If not, then remove that chr from the 'position'
position<-c(1:chr1,1:chr2,1:chr3,1:chr4,1:chr5,1:chr6,1:chr7)

#Change chr names: "1H" for "1"

CHR_names<-c(rep(1,chr1),rep(2,chr2),rep(3,chr3),rep(4,chr4),rep(5,chr5),rep(6,chr6),rep(7,chr7))

#Convert cM into Morgans
MORGANS<-GENOTYPE_READY$cm/100
SNP_file<-cbind(as.data.frame(GENOTYPE_READY[,1]),as.data.frame(CHR_names), as.data.frame(MORGANS),as.data.frame(position), as.data.frame(GENOTYPE_READY[,6]), as.data.frame(Major_Allele_ref))
colnames(SNP_file)<-c("SNP_name","Chromosome","cM","Position","Reference_mino","Reference_major")

#write.table(SNP_file,"~/Documents/Anita/THESIS/Data/for_smartPCA/Infiles_land_2/land_4527.snp",quote=F,row.names=F,col.names=F,sep="\t")

## 2. Create Sample.ind
Samples_ind<-cbind(as.data.frame(colnames(genotypes_or)), as.data.frame(rep('U',dim(genotypes_or)[2])), as.data.frame(rep(1,dim(genotypes_or)[2])))
colnames(Samples_ind)<-c("SampleID","Gender","pop_group")

#write.table(Samples_ind,"~/Documents/Anita/THESIS/Data/for_smartPCA/Infiles_land_2/land_4527.ind",quote=F,row.names=F,col.names=F,sep="\t")

## 3. Create Genotype.eigenstratgeno
COUNT_ALLELE<-function(dat){
	GENOTYPES<-dat[-c(1:6)]
	MAF<-dat[6]
	
	AA<-which(GENOTYPES == 'AA')
	BB<-which(GENOTYPES == 'BB')
	AB<-which(GENOTYPES == 'AB')
	Missing<- which(is.na(GENOTYPES) == 'TRUE')
	
	if(MAF == 'A'){
		GENOTYPES[AA]<-'2'
		GENOTYPES[BB]<-'0'
	}
	if(MAF == 'B'){
		GENOTYPES[AA]<-'0'
		GENOTYPES[BB]<-'2'
	}
	
	#All AB are 1
	GENOTYPES[AB]<-1
	
	#All missing values are = 9
	GENOTYPES[Missing]<-'9'
	return(GENOTYPES)
}

RESULT<-as.data.frame(apply(GENOTYPE_READY,1, COUNT_ALLELE))

#Turn table to have SNPs in rows and samples in columns

Genotype_counts<-as.data.frame(t(RESULT))

row.names(Genotype_counts)<-GENOTYPE_READY$SNP
colnames(Genotype_counts)<-colnames(GENOTYPE_READY[,7:(dim(GENOTYPE_READY)[2])])

write.table(Genotype_counts,"~/Documents/Anita/THESIS/Data/for_smartPCA/Infiles_land_2/land_4527.eigenstratgeno",quote=F,row.names=F,col.names=F,sep="")

#Now use the files as input files for smartPCA in the command line in MSI. package part of Eigensatat.Using a par file containing the parameters used
#E.g. smartpca -p par_land_4527 >./logfile_land_4527

#Use Plot_smartPCA.R to process output