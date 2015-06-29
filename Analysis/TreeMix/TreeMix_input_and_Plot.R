#Call genotype files for landraces and wild. This files contain the SNPs shared between land and WBDC. Snps are in columns, and samples in rows
rm(list=ls())

genotype_land<-read.table("~/Documents/github/BarleyLandraces/Datasets/Land_6152_SNPs_AB.txt",header=T,row.names=1)
genotype_wbdc<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/Wild_1896SNPs_AB.txt",header=T,row.names=1)

#Identify SNPs shared in both datasets
Shared_SNP<-intersect(colnames(genotype_land),colnames(genotype_wbdc))

#Select from both datasets only the SNP shared and order them by SNP name
land_snp<-genotype_land[,(colnames(genotype_land) %in% Shared_SNP)]
land<-land_snp[,order(colnames(land_snp))]

wbdc_snp<-genotype_wbdc[,(colnames(genotype_wbdc) %in% Shared_SNP)]
wbdc_snp_or<-wbdc_snp[,order(colnames(wbdc_snp))]
#Replace '??' for NA

Replace_NA<-function(dat){
	dat[which(dat == '??')]<-NA
	return(dat)
}
wbdc<-as.data.frame(apply(wbdc_snp_or,2,Replace_NA))
# Call population assignment results.land=k4. wbdc=k6

land_pop<-read.table("~/Documents/Structure/land_assignment_k4_iSelect6135.txt",header=T)
head(land_pop)
wbdc_pop<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)
head(wbdc_pop)

land[1:10,1:5];dim(land)
wbdc[1:10,1:5];dim(wbdc)

#Separate data based on populations
# For landraces
Pop1_l<-subset(land_pop,land_pop$Populations == '1')
Pop2_l<-subset(land_pop,land_pop$Populations == '2')
Pop3_l<-subset(land_pop,land_pop$Populations == '3')
Pop4_l<-subset(land_pop,land_pop$Populations == '4')

Pop1_land<-land[(row.names(land) %in% Pop1_l$Samples),]
Pop2_land<-land[(row.names(land) %in% Pop2_l$Samples),]
Pop3_land<-land[(row.names(land) %in% Pop3_l$Samples),]
Pop4_land<-land[(row.names(land) %in% Pop4_l$Samples),]

#For WBDC
Pop1_w<-subset(wbdc_pop,wbdc_pop$pop6_hap == '1')
Pop2_w<-subset(wbdc_pop,wbdc_pop$pop6_hap == '2')
Pop3_w<-subset(wbdc_pop,wbdc_pop$pop6_hap == '3')
Pop4_w<-subset(wbdc_pop,wbdc_pop$pop6_hap == '4')
Pop5_w<-subset(wbdc_pop,wbdc_pop$pop6_hap == '5')
Pop6_w<-subset(wbdc_pop,wbdc_pop$pop6_hap == '6')

Pop1_wbdc<-wbdc[(row.names(wbdc) %in% Pop1_w$Sample),]
Pop2_wbdc<-wbdc[(row.names(wbdc) %in% Pop2_w$Sample),]
Pop3_wbdc<-wbdc[(row.names(wbdc) %in% Pop3_w$Sample),]
Pop4_wbdc<-wbdc[(row.names(wbdc) %in% Pop4_w$Sample),]
Pop5_wbdc<-wbdc[(row.names(wbdc) %in% Pop5_w$Sample),]
Pop6_wbdc<-wbdc[(row.names(wbdc) %in% Pop6_w$Sample),]

##The next option is when WBDC divided EAST WEST
Pop_w_east<-rbind(as.data.frame(Pop1_wbdc),as.data.frame(Pop4_wbdc),as.data.frame(Pop6_wbdc))
Pop_w_west<-rbind(as.data.frame(Pop2_wbdc),as.data.frame(Pop3_wbdc),as.data.frame(Pop5_wbdc))

###Count allele frequency at each SNP for each population
#Using diploid unphased data
#Snps are in columns and samples in rows. The output will have SNPs in ROWs.
POPULATIONS<-c("Pop1_wbdc", "Pop2_wbdc", "Pop3_wbdc", "Pop4_wbdc", "Pop5_wbdc", "Pop6_wbdc","Pop1_land", "Pop2_land", "Pop3_land", "Pop4_land")
CountPop <-matrix(NA,nrow = c((dim(land)[2])),ncol=length(POPULATIONS))
colnames(CountPop)<-c(POPULATIONS)
for (p in 1:length(POPULATIONS)){
	dat <- get(POPULATIONS[p])
	
	for (s in 1:(dim(dat)[2])){
		AA<-length(which(dat[,s] == 'AA'))
		BB<-length(which(dat[,s] == 'BB'))
		AB<-length(which(dat[,s] == 'AB'))
		CountA<-(2*AA)+AB
		CountB<-(2*BB)+AB
		Table_count<-cbind(CountA, CountB)
		CountPop[s,p]<-paste(Table_count,collapse=",")
	}
}
colnames(CountPop)<-c("W_N_Zagros","W_S_Levant","W_S_Desert","W_Caspian_Sea","W_N_levant","W_C_asia","Centra_European","Asian","Coastal_Mediterranean","East_Africa")
write.table(CountPop,"~/Documents/TreeMix/Input/TreeMix_input.gz",quote=F,row.names=F,col.names=T,sep="\t")

#Add SNP names


row.names(CountPop)<-colnames(land)
write.table(CountPop,"~/Documents/TreeMix/Input/TreeMix_input_SNPnames.gz",quote=F,row.names=F,col.names=T,sep="\t")


#RUN TREEMIX
#treemix -i TreeMix_input.gz -o my_output

#PLOT TREE
source("~/Scripts/TreeMix/src/plotting_funcs.R")
plot_tree("~/Documents/TreeMix/output/my_output")
