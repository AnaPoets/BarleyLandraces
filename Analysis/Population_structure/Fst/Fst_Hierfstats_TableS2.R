rm(list=ls())

library(hierfstat)

#Load the genetic assignment for k=4 landraces (Clumpp_out.R)
assig<-read.table("~/Documents/land_assignment_k4_iSelect6152.txt",header=T)

#Load input file for STRUCTURE samples in rows, markers=Columns (Structure_1652_input.R)
Structure_input<-read.table("~/Documents/Str_input_1652.txt",header=T,row.names=1)

Structure_input <-Structure_input[,-c(1,2)]

#Change -9 for NA
Change_9<-function(dat){
	dat[which(dat == "-9")]<-NA
	return(dat)
}
genotypes<-as.data.frame(apply(Structure_input,2, Change_9))
#Select markers
markers<-colnames(genotypes)

#Select samples per population

head(assig)

pop1<-subset(assig, assig$Populations == "1")
pop2<-subset(assig, assig$Populations == "2")
pop3<-subset(assig, assig$Populations == "3")
pop4<-subset(assig, assig$Populations == "4")


### the next tree rows will arrange the files for Focal Fst analysis
#Central European vs. All
rbind(pop4,pop2,pop3)->all_pops
genotypes[(row.names(genotypes) %in% all_pops$Samples),]->TOTAL
genotypes[(row.names(genotypes) %in% pop1$Samples),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
	Fst_1_vsAll<-Fst

#Asian vs. All
rbind(pop1,pop4,pop3)->all_pops
genotypes[(row.names(genotypes) %in% all_pops$Samples),]->TOTAL
genotypes[(row.names(genotypes) %in% pop2$Samples),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
	Fst_2_vsAll<-Fst

#Coastal Mediterranean vs. All
rbind(pop1,pop2,pop4)->all_pops
genotypes[(row.names(genotypes) %in% all_pops$Samples),]->TOTAL
genotypes[(row.names(genotypes) %in% pop3$Samples),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
	Fst_3_vsAll<-Fst

#East Africa vs. all
rbind(pop1,pop2,pop3)->all_pops
genotypes[(row.names(genotypes) %in% all_pops$Samples),]->TOTAL
genotypes[(row.names(genotypes) %in% pop4$Samples),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
	Fst_4_vsAll<-Fst


#Make table S2

TableS2<-matrix(ncol=2,nrow=4)
colnames(TableS2)<-c("Median Fst","Maximum Fst")
row.names(TableS2)<-c("Central European","Asian","Coastal Mediterranean","East African")

for ( i in 1:4){
	TableS2[i,1]<-median(get(paste("Fst_",i,"_vsAll",sep="")),na.rm=T)
	TableS2[i,2]<-max(get(paste("Fst_",i,"_vsAll",sep="")),na.rm=T)
}

write.table(TableS2,"~/Documents/Additionalfile4_TableS2.xls",quote=F,row.names=T,col.names=T,sep="\t")