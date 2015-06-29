#Create input file for STRUCTURE, Using all the 803 landrace accessions with all 6,152SNPs

rm(list=ls())

land<-read.table("~/Documents/github/BarleyLandraces/Datasets/Land_6152_SNPs_AB.txt")

#Replace AA=1, BB=2, AB=-9,NA=-9

ReplaceGenotype<-function(dat){
	dat[which(dat == 'AA')] <-1
	dat[which(dat == 'BB')]<-2
	dat[which(dat == 'AB')] <-"-9"
	dat[which(dat == 'BA')] <-"-9"
	dat[which(is.na(dat))]<-"-9"
	return(dat)
}

genotypes_changed<-as.data.frame(apply(land,2, ReplaceGenotype))

#Add two columns one with 1s and one with 0s. The 0 indicates that there are not pre-assign individuals.

Str_input<-cbind(rep(7,(dim(genotypes_changed)[1])),rep(0,(dim(genotypes_changed)[1])), genotypes_changed)

write.table(Str_input,"~/Desktop/Str_input_1652.txt",quote=F,row.names=T,col.names=F,sep="\t")