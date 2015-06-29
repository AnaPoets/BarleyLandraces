#After runing IBS_input_PEDfiles_1896SNPs.R
#THIS CODE WILL:
#1. Remove the SNPs at the end of each chromosome that can not be devided in 30 (30SNP windows)
#2. Create new map files with only the remaining SNPs
#3.	Split each chromosome in 30SNP windows for both the new map and the ped files.

rm(list=ls())

#1. .map files created in IBS_input_PEDfiles_1896SNPs and moved to the IBS directory
chr1<-read.table("~/Documents/Input/ByChromosome/chr1.map")
chr2<-read.table("~/Documents/Input/ByChromosome/chr2.map")
chr3<-read.table("~/Documents/Input/ByChromosome/chr3.map")
chr4<-read.table("~/Documents/Input/ByChromosome/chr4.map")
chr5<-read.table("~/Documents/Input/ByChromosome/chr5.map")
chr6<-read.table("~/Documents/Input/ByChromosome/chr6.map")
chr7<-read.table("~/Documents/Input/ByChromosome/chr7.map")


for (i in 1:7){
	#Set here the window's size you want to use.
	SNPcount <-30
	CHR<-get(paste("chr",i,sep=""))
	CHR_SNP<-floor(dim(CHR)[1]/SNPcount)* SNPcount
	#Write .map files with just the number of SNPs that are devided by 30.
	Input<-as.data.frame(CHR[1:CHR_SNP,])
	Output<-paste("~/Documents/Anita/THESIS/Data/for_IBS/Barley_1896/Input/PLINK/chr",i,"_30s.map",sep="")
	write.table(Input,Output,,quote=F,row.names=F,col.names=F,sep="\t")
	}

###2. Devide the map and .ped files in segments of 30SNPs


for (c in 1:7){

CHR<-c
Chr_map<-read.table(paste("~/Documents/Anita/THESIS/Data/for_IBS/Barley_1896/Input/PLINK/chr",CHR,"_30s.map",sep=""),header=F)
W_L_ped<-read.table(paste("~/Documents/Anita/THESIS/Data/for_IBS/Barley_1896/Input/ByChromosome/Wild_land_chr",CHR,".ped",sep=""),header=F)

#Separate sample names from ped file

SAMPLES<-W_L_ped[,1:6]
W_L_ped<-W_L_ped[,-c(1:6)]

#Select the same number of markers that the new 30s.map has. Remember that the .ped has two columns per marker
#since we removed the end of each chromosome, we can select only the lenght that reminds times 2, for two alleles per marker
Trim_ped<-W_L_ped[,1:(dim(Chr_map)[1]*2)] 
dim(Trim_ped)

array_ped_30s<-seq(1, (dim(Trim_ped)[2]),(SNPcount*2))
array_map_30s<-seq(1, (dim(Chr_map)[1]), SNPcount)

	#Separate map in 30SNP map files and .ped files.
	
	for (i in 1:length(array_map_30s)){
		Map30<-Chr_map[(array_map_30s[i]):(array_map_30s[i]+(SNPcount-1)),]
		Ped30<-Trim_ped[,(array_ped_30s[i]): (array_ped_30s[i]+((SNPcount*2)-1))]
		Ped30_ready<-cbind(as.data.frame(SAMPLES),as.data.frame(Ped30))
		OUTPUT_map<-paste("~/Documents/PLINK/output/Sep_30SNP/Chr",c,"_",i,".map",sep="")
		OUTPUT_ped<-paste("~/Documents/PLINK/output/Sep_30SNP/Chr",c,"_",i,".ped",sep="")
		write.table(Map30,OUTPUT_map,quote=F,row.names=F,col.names=F,sep="\t")
		write.table(Ped30_ready,OUTPUT_ped,quote=F,row.names=F,col.names=F,sep="\t")
		
		#Make a list of the Starting and End Snp for each chromosomal segment of 30SNPs
		Segments_positions<-Chr_map[c(array_map_30s[i]),]
		Segments_info<-cbind(paste("Chr",c,"_",i,sep=""),c("start"), Chr_map[array_map_30s[i],])
		Segments_info2<-cbind(paste("Chr",c,"_",i,sep=""),c("end"), Chr_map[c(array_map_30s[i]+(SNPcount-1)),])
		SEGMENTS_INFORMATION<-cbind(Segments_info, Segments_info2)
		
		write.table(SEGMENTS_INFORMATION,"~/Documents/PLINK/Sep_30SNP/Segments_info.txt",append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")
		
		write.table(Map30,"~/Documents/PLINK/Sep_30SNP/map_cm_30snpALL.map",append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")
		
		}
	
}

#Run output in plink: Step3_Plink_run.sh
#in the terminal make a list of all the file names without extension 
#find . -type f -name "*ped" | sed 's/\.\///g' | sed 's/.ped//g' >./List_gral.txt