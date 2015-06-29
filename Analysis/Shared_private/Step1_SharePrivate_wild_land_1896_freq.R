##Using Private SNPs to each wild population
##This code looks to identify which individual landrace carry which private allele from the wild.
##Use diploid unphased data. Any heterozygote site would be counted on favor or carring the minor allele from the wild


rm(list=ls())


genotypes_land<-read.table("~/Documents/github/BarleyLandraces/Datasets/Land_1896SNPs_AB.txt",header=T,row.names=1)

private_alleles<-read.table("~/Documents/github/BarleyLandraces/Datasets/Private_allele_ALLpops.txt",header=T)

head(private_alleles)

#Make Minor Allele state diploid
MinorAll_diploid<-function(dat){
	dat[dat == 'B']<-'BB'
	dat[dat == 'A']<-'AA'
	return(dat)
}
private<-as.data.frame(apply(private_alleles,2, MinorAll_diploid))
#order private alleles by name
private_or<-private[order(private[,1]),]

#Select land SNPs that are private in wild
land2<-land[,(colnames(land) %in% private_or[,1])]
dim(land2)
land2[1:10,1:10]

#turn the table to add the minor allele 
as.data.frame(t(land2))->t_land2
t_land2[1:10,1:10]
dim(t_land2)
#Order land shared alleles by SNP name
land2_or<-t_land2[order(row.names(t_land2)),]
dim(land2_or)
head(land2_or[1:10,1:10])
sample_name<-colnames(land2_or)

#Add minor_allele to land2_or
#Check that the SNPs are in the same order
if (identical(private_or[,"SNP"],as.data.frame(row.names(land2_or))[,"row.names(land2_or)"])) {
minor_allele<-private_or[,2]
land_minor<-cbind(as.data.frame(minor_allele), land2_or)
land_minor[1:10,1:10]
} else print ("ERROR! SNP names don't match")

#Put minor allele at the top of the file. Columns are SNPs , rows are individuals
t_land_minor<-t(land_minor)
t_land_minor<-as.data.frame(t_land_minor)
t_land_minor[1:10,1:10]

#Separate private by wild pop
N_mesopotamia<-subset(private_or, private_or $Population == 'N_Mesopotamia')
S_Levant<-subset(private_or, private_or $Population == 'S_Levant')
S_desert<-subset(private_or, private_or $Population == 'S_Desert')
N_levant<-subset(private_or, private_or $Population == 'N_Levant')
C_Asia<-subset(private_or, private_or $Population == 'C_Asia')



#Separate the landraces' markers according to SNPs private to each wild population
land_N_meso<-t_land_minor[,(colnames(t_land_minor) %in% N_mesopotamia[,1])]
land_S_Levant <-t_land_minor[,(colnames(t_land_minor) %in% S_Levant[,1])]
land_S_desert <-t_land_minor[,(colnames(t_land_minor) %in% S_desert[,1])]
land_N_levant <-t_land_minor[,(colnames(t_land_minor) %in% N_levant[,1])]
land_C_Asia<-t_land_minor[,(colnames(t_land_minor) %in% C_Asia[,1])]


#Frequency of N_mesopotamia alleles in landraces. Identify which individuals carry the minor allele at each of these private alleles.

results_N_meso<-matrix(NA,dim(land_N_meso)[1],dim(land_N_meso)[2])
#FOR each individual i , look at each SNP j .If private allele present then 1, otherwise 0
for(i in 2:(dim(land_N_meso)[1])){
	for (j in 1:(dim(land_N_meso)[2])) 
if(land_N_meso[i,j] == land_N_meso[1,j]) {'1'} -> results_N_meso[i,j]
	else if (land_N_meso[i,j] == "AB") {'11'} -> results_N_meso[i,j]
	else {"0"} -> results_N_meso[i,j]
}
land_N_meso_freq<-as.data.frame(results_N_meso[-1,])
colnames(land_N_meso_freq)<-colnames(land_N_meso)
row.names(land_N_meso_freq)<-row.names(land)

#Frequency of S_levant alleles
results_S_levant<-matrix(NA,dim(land_S_Levant)[1],dim(land_S_Levant)[2])
#FOR each individual i , look at each SNP j
for(i in 2:(dim(land_S_Levant)[1])){
	for (j in 1:(dim(land_S_Levant)[2]) )
if(land_S_Levant[i,j] == land_S_Levant[1,j]) {'2'} -> results_S_levant[i,j]
	else if (land_S_Levant[i,j] == "AB") {'22'} -> results_S_levant[i,j]
	else {"0"}-> results_S_levant[i,j]
}
land_S_levant_freq<-as.data.frame(results_S_levant[-1,])
colnames(land_S_levant_freq)<-colnames(land_S_Levant)
row.names(land_S_levant_freq)<-row.names(land)

#Frequency of S_desert alleles
results_S_desert<-matrix(NA,dim(land_S_desert)[1],dim(land_S_desert)[2])
#FOR each individual i , look at each SNP j
for(i in 2:(dim(land_S_desert)[1])){
	for (j in 1:(dim(land_S_desert)[2]) )
if(land_S_desert[i,j] == land_S_desert[1,j]) {'3'} -> results_S_desert[i,j]
	else if (land_S_desert[i,j] == "AB") {'33'} -> results_S_desert[i,j]
	else {"0"}-> results_S_desert[i,j]
}
land_S_desert_freq<-as.data.frame(results_S_desert[-1,])
colnames(land_S_desert_freq)<-colnames(land_S_desert)
row.names(land_S_desert_freq)<-row.names(land)

#Frequency of N_levant alleles
results_N_levant<-matrix(NA,dim(land_N_levant)[1],dim(land_N_levant)[2])
#FOR each individual i , look at each SNP j
for(i in 2:(dim(land_N_levant)[1])){
	for (j in 1:(dim(land_N_levant)[2]) )
if(land_N_levant[i,j] == land_N_levant[1,j]) {'5'} -> results_N_levant[i,j]
	else if (land_N_levant[i,j] == "AB") {'55'} -> results_N_levant[i,j]
	else {"0"}-> results_N_levant[i,j]
}
land_N_levant_freq<-as.data.frame(results_N_levant[-1,])
colnames(land_N_levant_freq)<-colnames(land_N_levant)
row.names(land_N_levant_freq)<-row.names(land)


#Frequency of C_Asia alleles
results_C_Asia<-matrix(NA,dim(land_C_Asia)[1],dim(land_C_Asia)[2])
#FOR each individual i , look at each SNP j
for(i in 2:(dim(land_C_Asia)[1])){
	for (j in 1:(dim(land_C_Asia)[2]) )
if(land_C_Asia[i,j] == land_C_Asia[1,j]) {'6'} -> results_C_Asia[i,j]
	else if (land_C_Asia[i,j] == "AB") {'66'} -> results_C_Asia[i,j]	
	else {"0"}-> results_C_Asia[i,j]
}
land_C_Asia_freq<-as.data.frame(results_C_Asia[-1,])

colnames(land_C_Asia_freq)<-colnames(land_C_Asia)
row.names(land_C_Asia_freq)<-row.names(land)

#Combine the results

AncestryPrivate_present<-cbind(as.data.frame(land_N_meso_freq),as.data.frame(land_S_levant_freq),as.data.frame(land_S_desert_freq),as.data.frame(land_N_levant_freq),as.data.frame(land_C_Asia_freq))

write.table(AncestryPrivate_present,"~/Documents/SharedPoly/Land_wbdc_1896/Frequency_private_alleles/Private_PresentInLand.txt",quote=F,row.names=T,col.names=T,sep="\t")

#Continue with SharedPrivate_alleleFreq_1896SNPs_diploid.R