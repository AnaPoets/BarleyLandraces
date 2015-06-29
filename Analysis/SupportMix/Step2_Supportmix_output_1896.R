rm(list=ls())

#Run the first lines until write.table for each chromosome
#Create a subdirectory to keep the output from chr_good

dir.create(paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob",sep=""))

for (i in 1:7) {
#Import population assigment from SupportMix
infile_pops<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/outSupportMix_chr",i,'.tped',sep="")
pops<-read.table(infile_pops)

#Import the probability of assignment
infile_probs<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/outSupportMix_chr",i,'.Probs.tped',sep="")
probs<-read.table(infile_probs)
pops[,1:10]
probs[,1:10]

#Matrix to hold the ancestry assignments, after puting "6" for sites with less than <95% prob of assignment
GOOD<-matrix(ncol=(dim(pops)[2]),nrow=(dim(pops)[1]))

#by columns= samples without the SNP information. For each ancestry assignment at each individual haplotype, look the prob of assignment at each SNP window.
for(j in 5:(dim(pops)[2])) {
	#for each SNP window at each sample
	for(s in 1:(dim(pops)[1])){
		if ((probs[s,j]) < 0.95) GOOD[s,j] <- "6" else GOOD[s,j] <- pops[s,j]
	}
}

#Remove the first 4 columns of SNP information that are empty and add the right information, Samples are in COLUMNS and SNP windows in ROWS

Good<-GOOD[,-c(1:4)]
Good<-as.data.frame(t(cbind(pops[,1:4],Good)))

#in the created a folder called Pop_prob
output<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr",i,"_good",sep="")
write.table(Good,output,sep="\t",quote=F,row.names=F,col.names=F)
}

#Call the corrected files
chr1_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr1_good")

chr2_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr2_good")

chr3_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr3_good")

chr4_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr4_good")

chr5_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr5_good")

chr6_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr6_good")

chr7_good<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/pop_prob/chr7_good")

#Import samples names

samples_n<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/tpedData/land_admixed_k5.tfam.gz",header=F)
Index_odd<-seq(1,2*(dim(samples_n)[1]),2)
head(Index_odd)
names_odd<-as.data.frame(cbind(as.data.frame(samples_n[,1]),Index_odd))
head(names_odd)
colnames(names_odd)<-c("Samples","Index")

samples_even<-as.data.frame(paste(samples_n[,1],"_2",sep=""))
Index_even<-seq(2,2*(dim(samples_n)[1]),2)
names_even<-as.data.frame(cbind(as.data.frame(samples_even[,1]), Index_even))
head(names_even)
colnames(names_even)<-c("Samples","Index")

col_names<-as.data.frame(rbind(names_even,names_odd))
head(col_names)
col_names_or<-col_names[order(col_names$Index),]
head(col_names_or)
col_names_or_samples<-as.data.frame(col_names_or[,1])
colnames(col_names_or_samples)<-c("Info")
Support_results_all<-cbind(as.data.frame(chr1_good),as.data.frame(chr2_good),as.data.frame(chr3_good),as.data.frame(chr4_good),as.data.frame(chr5_good),as.data.frame(chr6_good),as.data.frame(chr7_good))
dim(Support_results_all)
#Turn table to have SNPwindows in ROWs and Samples in columns
Support_results<-as.data.frame(t(Support_results_all))

dim(Support_results)

SNP_infoNAMES<-as.data.frame(c("Chro","SNP","Cumulative","Extra"))
colnames(SNP_infoNAMES)<-c("Info")

col_names_all<-rbind(SNP_infoNAMES, col_names_or_samples)

colnames(Support_results)<-col_names_all[,1]
Support_results[1:10,1:10]
#Add an index column to keep the chromosomes and SNPs in order

Index<-c(1:dim(Support_results)[1])

support_indexed<-cbind(as.data.frame(Index),as.data.frame(Support_results))
head(support_indexed[1:10,1:10])

###================================================
#call iSelect genetic map
genetic_map<-read.table("~/Documents/github/BarleyLandraces/Datasets/GeneticMap_iSelect_9k.txt",header=T)
head(genetic_map)

#select SNPs in genetic map that were used as markers for supportMix

genetic_map_support<-genetic_map[(genetic_map[,1] %in% support_indexed$SNP),]
head(genetic_map_support)

#Paste gmap and support results

gmap_and_support<-cbind(as.data.frame(genetic_map_support),as.data.frame(support_indexed))
gmap_and_support[1:10,1:15]
###calculate the interval between the SNPs in the table gmap_and_support and add this column to the file in the far left 'Interval'

output_gmap<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/gmap_and_suport.xls",sep="")
write.table(gmap_and_support, output_gmap,quote=F,sep="\t",col.names=T,row.names=F)


#change pop1(0)=1,pop2(1)=2, pop3(2)=3,pop5(3)=5,pop6(4)=6, Unknown(6)=7
SNP_info<-gmap_and_support[,c(1:12)]
Indiv_assig<-gmap_and_support[,-c(1:12)]

change_assig<-function(dat){
	dat[dat == '0'] <-'A'
	dat[dat == '1'] <-'B'
	dat[dat == '2'] <-'C'
	dat[dat == '3'] <-'D'
	dat[dat == '4'] <-'E' 
	dat[dat == '6'] <-'F'
	return(dat)
}
assign_letters<-as.data.frame(apply(Indiv_assig,2,change_assig))

change_assig_again<-function(dat){
	dat[dat == 'A'] <-'1'
	dat[dat == 'B'] <-'2'
	dat[dat == 'C'] <-'3'
	dat[dat == 'D'] <-'5'
	dat[dat == 'E'] <-'6' 
	dat[dat == 'F'] <-'7'
	return(dat)
}
assign_col<-as.data.frame(apply(assign_letters,2,change_assig_again))

#Add the Interval column
Cumulative<-SNP_info$Cumulative
Cumulative[length(Cumulative) +1] <- '1112.7'
Cumulative<-as.data.frame(Cumulative)
Interval<-(rep(NA,(dim(Cumulative)[1]-1)))
for (i in 1:dim(Cumulative)[1]-1) {Interval[i]<-as.numeric(as.matrix(Cumulative[(i+1),])) - as.numeric(as.matrix(Cumulative[i,]))}

#Put together Interval, SNP infor and the assignation with the numbers changed to have right colors
gmap_and_suport_complete<-cbind(as.data.frame(Interval),as.data.frame(SNP_info),as.data.frame(assign_col))

write.table(gmap_and_suport_complete,"~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/gmap_and_suport_complete.xls",quote=F,sep="\t",col.names=T,row.names=F)


##NOW CONTINUE WITH CODE "Step3_SupportMix_calculate_proportions.R"