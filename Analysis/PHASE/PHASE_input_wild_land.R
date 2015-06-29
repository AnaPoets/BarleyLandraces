rm(list=ls())

land<-read.table("~/Documents/github/BarleyLandraces/Datasets/Land_6152_SNPs_AB.txt",header=T,row.names=1)
wild<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/Wild_1896SNPs_AB.txt",header=T,row.names=1)

#Remove Caspian Sea wild samples
wild_pops<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)
caspian<-subset(wild_pops, wild_pops$pop6_hap == '4')
wild<-wild[!(row.names(wild) %in% caspian$Sample),]

land[1:10,1:10]
wild[1:10,1:10]

SNPs_common<-intersect(colnames(land),colnames(wild))

#Identify how many of those SNPs have genetic position

genmap<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)

head(genmap)

SNPs_common_map<-intersect(SNPs_common,genmap$SNP)
length(SNPs_common_map)


#Select SNPs from land and wild

land_snp<-land[,(colnames(land) %in% SNPs_common_map)]

wild_snp<-wild[,(colnames(wild) %in% SNPs_common_map)]

#Order SNPs according to genetic map

genmap1896<-genmap[(genmap$SNP %in% SNPs_common),]

list_snp_shared_land_wild<-genmap1896[,c(2,3,4)]

genmap1896_or<-genmap1896[order(genmap1896$SNP),]


#rotate land and wild and order them by SNP name
land_t<-t(land_snp)
wild_t<-t(wild_snp)

land_t_or<-land_t[order(row.names(land_t)),]
wild_t_or<-wild_t[order(row.names(wild_t)),]

#Add SNPs from genetic map

land_genmap<-cbind(as.data.frame(genmap1896_or),as.data.frame(land_t_or))
wild_genmap<-cbind(as.data.frame(genmap1896_or),as.data.frame(wild_t_or))

#Sort by Index

land_genmap_or<-land_genmap[order(land_genmap$index),]

wild_genmap_or<-wild_genmap[order(wild_genmap$index),]

land_genmap_or[1:10,1:10]

#Remove extra columns and turn tables to have SNPs in columnas and samples in rows

land_genmap_or_ready<-land_genmap_or[,-c(2:7)]
land_genmap_or_ready_t<-t(land_genmap_or_ready)

wild_genmap_or_ready<-wild_genmap_or[,-c(2:7)]
wild_genmap_or_ready_t<-t(wild_genmap_or_ready)

#Make the first row the column headers
colnames(land_genmap_or_ready_t)<-land_genmap_or_ready_t[1,]
land_genmap_or_ready_t<-as.data.frame(land_genmap_or_ready_t[-1,])
land_genmap_or_ready_t[1:10,1:10]

colnames(wild_genmap_or_ready_t)<-wild_genmap_or_ready_t[1,]
wild_genmap_or_ready_t <-as.data.frame(wild_genmap_or_ready_t[-1,])
wild_genmap_or_ready_t[1:10,1:10]



######################################
##Separate Genotypes AA in A A, BB in B B, AB in A B, NA in ? ?
#####Select landraces or wild to work on land_genmap_or_ready_t or wild_genmap_or_ready_t

POPULATIONS<-c("land_genmap_or_ready_t","wild_genmap_or_ready_t")
Populations_abbreviation<-c("land","wbdc")
for (p in 1:lenght(POPULATIONS)){
DATA<-get(POPULATIONS[p])

SNPnames<-colnames(DATA)
Sep_genotypes<-function(dat){
	dat [dat == 'AA'] <-c('A	A')
	dat [dat == 'BB'] <-c('B	B')
	dat [dat == 'AB'] <-c('A	B')
	dat [is.na(dat)] <-c('?	?') 
	return (dat)
}
DOUBLED_GENOTYPES<-t(as.data.frame(apply(DATA,1,Sep_genotypes)))
write.table(DOUBLED_GENOTYPES,"~/Desktop/temp.txt",quote=F,row.names=T,col.names=F,sep="\t")

#Import file again without header so we can split the genotypes and add SNPs names

Genotypes<-read.table("~/Desktop/temp.txt",header=F,row.names=1)
SAMPLE<-row.names(Genotypes)
#Create a string of Sample names that had and extra _2, so we can merge the even and odd genotypes and get only one file 

SAMPLE_DOBLE_LAND <-as.data.frame(paste(SAMPLE,"_2",sep=""))

#Separate genotypes in haplotypes. Add SNPs and Samples_2 names
ODD_columns<-Genotypes[,seq(1,ncol(Genotypes),by=2)]
colnames(ODD_columns) <-SNPnames
EVEN_columns<-Genotypes[,seq(2,ncol(Genotypes),by=2)]
row.names(EVEN_columns) <-SAMPLE_DOBLE_LAND [,1]
colnames(EVEN_columns)<-SNPnames

#Combine both files

DOUBLED_DATA<-rbind(as.data.frame(ODD_columns), as.data.frame(EVEN_columns))

#Sort by sample name

DOUBLED_DATA_or <-DOUBLED_DATA[order(row.names(DOUBLED_DATA)),]

#Devide samples by chromosome
Chr1<-subset(genmap1896, genmap1896$Chro == '1H')
Chr2<-subset(genmap1896, genmap1896$Chro == '2H')
Chr3<-subset(genmap1896, genmap1896$Chro == '3H')
Chr4<-subset(genmap1896, genmap1896$Chro == '4H')
Chr5<-subset(genmap1896, genmap1896$Chro == '5H')
Chr6<-subset(genmap1896, genmap1896$Chro == '6H')
Chr7<-subset(genmap1896, genmap1896$Chro == '7H')

DOUBLED_DATA_or_chr1<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr1$SNP)]
DOUBLED_DATA_or_chr2<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr2$SNP)]
DOUBLED_DATA_or_chr3<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr3$SNP)]
DOUBLED_DATA_or_chr4<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr4$SNP)]
DOUBLED_DATA_or_chr5<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr5$SNP)]
DOUBLED_DATA_or_chr6<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr6$SNP)]
DOUBLED_DATA_or_chr7<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr7$SNP)]


for (i in 1:7) {
	INPUT<-get(paste("DOUBLED_DATA_or_chr",i,sep=""))
	OUTPUT<-paste("~/Documents/",Populations_abbreviation[p],"_prephase_chr_",i,".txt",sep="")

write.table(INPUT,OUTPUT,quote=F,row.names=T,col.names=T,sep="\t")
}
}
##Now run New_make_fasta.pl

