rm(list=ls())
library(LDheatmap)
library(genetics)

land<-read.table("~/Documents/github/BarleyLandraces/Datasets/Land_6152_SNPs_AB.txt",header=T,row.names=1)

#Import genetic map to separate SNP by Chromosome

genmap<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)

chr1<-subset(genmap,genmap$Chro == '1H')
chr2<-subset(genmap,genmap$Chro == '2H')
chr3<-subset(genmap,genmap$Chro == '3H')
chr4<-subset(genmap,genmap$Chro == '4H')
chr5<-subset(genmap,genmap$Chro == '5H')
chr6<-subset(genmap,genmap$Chro == '6H')
chr7<-subset(genmap,genmap$Chro == '7H')

#Separate Genotype by Chromosomes
geno_chr1<-land[,(colnames(land) %in% chr1$SNP)]
geno_chr2<-land[,(colnames(land) %in% chr2$SNP)]
geno_chr3<-land[,(colnames(land) %in% chr3$SNP)]
geno_chr4<-land[,(colnames(land) %in% chr4$SNP)]
geno_chr5<-land[,(colnames(land) %in% chr5$SNP)]
geno_chr6<-land[,(colnames(land) %in% chr6$SNP)]
geno_chr7<-land[,(colnames(land) %in% chr7$SNP)]

# Remove SNPs with Minor Allele Frequency (MAF) <0.05

low_freq <- function(dat) {
#try_size <- length(na.omit(try))
AA <- dat [dat == "AA"]
AA <- length(AA)
BB <- dat[dat == "BB"]
BB <- length(BB)
smaller <- min(c(AA,BB))

smaller/(sum(AA + BB)) ->SMALLER
return(SMALLER)
}

#For each Chromosome remove SNPs with MAF<0.05

for (chr in 1:7){
	CHR<-get(paste("geno_chr",chr,sep=""))
	low_freq_rm <- as.vector(which(apply(CHR,2,low_freq) < 0.05))
	if (length(low_freq_rm) > 0){
	CHR_new <-CHR[,c(-low_freq_rm)]	
	} else (CHR ->CHR_new)
	assign(paste("land_chr",chr,sep=""),CHR_new)
}

# transform genotypes separating in two haplotypes

genos1 <- makeGenotypes(land_chr1,sep='')
genos2 <- makeGenotypes(land_chr2,sep='')
genos3 <- makeGenotypes(land_chr3,sep='')
genos4 <- makeGenotypes(land_chr4,sep='')
genos5 <- makeGenotypes(land_chr5,sep='')
genos6 <- makeGenotypes(land_chr6,sep='')
genos7 <- makeGenotypes(land_chr7,sep='')

for ( i in 1:7){
pdf(paste("~/Documents/plot/LD_",i,"H_r2.pdf",sep=""),width=10,height=10)
LDheatmap(paste("genos",i,sep=""), LDmeasure="r", distances='genetic', add.map=F, flip=T)
dev.off()
}

