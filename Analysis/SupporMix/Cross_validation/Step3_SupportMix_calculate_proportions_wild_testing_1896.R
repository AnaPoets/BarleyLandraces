
# The following code will calculate the proportion of each Wild population for each landrace population by each chromosomal segment.
rm(list=ls())

#Define from which runs to analyze
INITIATE<-c(1)
FINALIZE<-c(50)
WindowsSize<-c("75_SNPs_output/")
#Set working directory's name
for (D in  INITIATE:FINALIZE){
DirectoryName<-paste("20_wild_testing_rep",D,sep="")


#DirectoryName<-c("20_wild_testing_rep7")
#r=7
####
#Create a subdirectory to keep the output from chr_good
dir.create(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob",sep=""))

for (i in 1:7) {
infile_pops<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/outSupportMix_chr",i,'.tped',sep="")
pops<-read.table(infile_pops)

#probs<-read.table("~/Dropbox/ANITA_301213/THESIS/Data/for_SupporMix/Supportmix_barley1846_ACTG_TO_AB/Wild_k5/100_windows/chr7/outSupportMix_chr7.Probs.tped")
#probs<-read.table("~/Documents/SupporMix/SupportMix_wild_testing_1896/75_SNPs/wild_run1_random/Analysis/outSupportMix_chr5.Probs.tped")

infile_probs<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/outSupportMix_chr",i,'.Probs.tped',sep="")
probs<-read.table(infile_probs)
pops[,1:10]
probs[,1:10]
GOOD<-matrix(ncol=(dim(pops)[2]),nrow=(dim(pops)[1]))

#by columns= samples without the SNP information
for(j in 5:(dim(pops)[2])) {
	#for each SNP window at each sample
	for(s in 1:(dim(pops)[1])){
		if ((probs[s,j]) < 0.95) GOOD[s,j] <- "6" else GOOD[s,j] <- pops[s,j]
	}
}

#Remove the first 4 columns of SNP information that are empty and add the right information, and turn table to have Samples in ROWS and SNP windows in columns
Good<-GOOD[,-c(1:4)]
Good<-cbind(pops[,1:4],Good)
chr_good<-t(Good)
#create a folder called pop_prob

output<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr",i,"_good",sep="")
write.table(chr_good,output,sep="\t",quote=F,row.names=F,col.names=F)
}

#Call the corrected files
chr1_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr1_good",sep=""))
chr2_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr2_good",sep=""))
chr3_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr3_good",sep=""))
chr4_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr4_good",sep=""))
chr5_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr5_good",sep=""))
chr6_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr6_good",sep=""))
chr7_good<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/pop_prob/chr7_good",sep=""))


#Import samples names

samples_n<-read.table(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Header_admixed_wbdc_testing.tped.gz", sep=""),header=T)
colnames(samples_n)->col_names

col_names_or_samples<-as.data.frame(col_names[-c(1:4)])

colnames(col_names_or_samples)<-c("Info")




Support_results_all <-cbind(as.data.frame(chr1_good),as.data.frame(chr2_good),as.data.frame(chr3_good),as.data.frame(chr4_good),as.data.frame(chr5_good),as.data.frame(chr6_good),as.data.frame(chr7_good))
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

#============
#call iSelect genetic map
genetic_map<-read.table("~/Dropbox/ANITA_301213/THESIS/Data/GeneticMap_iSelect_9k_012213.txt",header=T)
head(genetic_map)

#select SNPs in genetic map that were used as markers for supportMix

genetic_map_support<-genetic_map[(genetic_map[,1] %in% support_indexed$SNP),]
head(genetic_map_support)

#Paste gmap and support results

gmap_and_support<-cbind(as.data.frame(genetic_map_support),as.data.frame(support_indexed))
gmap_and_support[1:10,1:15]
###calculate the interval between the SNPs in the table gmap_and_support and add this column to the file in the far left 'Interval'

output_gmap<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/gmap_and_suport.xls",sep="")
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


output_complete<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/gmap_and_suport_complete.xls",sep="")
write.table(gmap_and_suport_complete, output_complete,sep="\t",quote=F,row.names=F,col.names=T)

##########Calculate proportions.R:

supmix_in<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/gmap_and_suport_complete.xls",sep="")
supmix<-read.table(supmix_in,header=T)
dim(supmix)

#Select the chrom windows information (cM, snp name, etc)
chr<-supmix[,1:13]
supmix_all_samples<-supmix[,-c(1:13)]
#call the population assignement according to Structure k=4.

str_wild<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)
head(str_wild)

#make another str file adding the '_2' to each sample name. Then combine the files to later separtate all the wild samples

str_wild_2<-cbind(paste(str_wild[,1],'_2',sep=""),as.data.frame(str_wild[,2:4]))
colnames(str_wild_2)<-c("Sample","pop6_hap" ,"latitude" ,"longitude")

assi_str_all<-rbind(as.data.frame(str_wild),as.data.frame(str_wild_2))
tail(assi_str_all)


#Create a subdirectory to keep the output from proportions
dir.create(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/Pop_summary",sep=""))

# Select samples names for each population

pop1<-subset(assi_str_all, assi_str_all $pop6_hap == "1")
pop2<-subset(assi_str_all, assi_str_all $pop6_hap == "2")
pop3<-subset(assi_str_all, assi_str_all $pop6_hap == "3")
pop4<-subset(assi_str_all, assi_str_all $pop6_hap == "4")
pop5<-subset(assi_str_all, assi_str_all $pop6_hap == "5")
pop6<-subset(assi_str_all, assi_str_all $pop6_hap == "6")
# Separate resutls from supportMix for each population

supmix[,(colnames(supmix) %in% pop1$Sample)]->supmix1
supmix[,(colnames(supmix) %in% pop2$Sample)]->supmix2
supmix[,(colnames(supmix) %in% pop3$Sample)]->supmix3
supmix[,(colnames(supmix) %in% pop5$Sample)]->supmix5
supmix[,(colnames(supmix) %in% pop6$Sample)]->supmix6

#Paste the chr information to each supmix population

cbind(as.data.frame(chr),as.data.frame(supmix1)) -> chr_supmix1
cbind(as.data.frame(chr),as.data.frame(supmix2)) -> chr_supmix2
cbind(as.data.frame(chr),as.data.frame(supmix3)) -> chr_supmix3
cbind(as.data.frame(chr),as.data.frame(supmix5)) -> chr_supmix5
cbind(as.data.frame(chr),as.data.frame(supmix6)) -> chr_supmix6

##CHANGE HERE THE POPULATION NUMBER THAT YOU WANT TO ANALYZE

#for each population (p)analyse the proportions

SUMMARY_w_unassigned<-matrix(ncol=6,nrow=6)
SUMMARY_all_assigned<-matrix(ncol=5,nrow=6)
for (p in c(1,2,3,5,6)) {

supmix_analyzed<-get(paste('supmix',p,sep=""))
t(supmix_analyzed)->t_supmix
dim(t_supmix)

#Calculate the number of times that each population in the Wild contributes to each segment (23) in the landraces populations.Accros populations what is the proportion at each segment.

matrix(0,dim(t_supmix)[[2]],7)->my_results

for (i in 1:dim(t_supmix)[2]){
	table(t_supmix[,i]) ->table_R
	as.data.frame(table_R)->my_table
	intersect(c(1,2,3,4,5,6,7),my_table[,1]) ->table_data
	as.data.frame(table_data)->items_value
	
	for (j in 1:dim(items_value)[1]) {
	items_value[j,1] ->position
	as.numeric(as.character(position))->position2
	     my_table[j,2] ->my_results[i,position2]}
	}
	
	my_results
	colnames(my_results)<-c("1","2","3","4","5","6","7")
	important_col<-my_results[,c(1,2,3,5,6,7)]

#Calculate the proportions of the contribution of each Wild (7 columns) to the "testing wild" population at each segment: 1)with Unassigned sites, 2)without unassigned
# 1) including unassined sites
proportions_w_unassignment <-matrix(0,dim(t_supmix)[[2]],6) ### activate this when calculating proportions including unknowns
for (i in 1:dim(t_supmix)[[2]]){
	important_col[i,]/apply(important_col,1,sum)[i] -> proportions_w_unassignment[i,]
}
dim(proportions_w_unassignment)
head(proportions_w_unassignment)
summary(proportions_w_unassignment)
#Add columns of information Interval, chr, Cumulative,index

cbind(as.data.frame(chr[,c(1,4,6,7)]), as.data.frame(proportions_w_unassignment)) ->Pop_summary_unassi
colnames(Pop_summary_unassi)<-c("Interval","Chro","index","Cumulative","red","blue","green","orange","purple","yellow")
Pop_summary_unassi ## remove yellow for just known things

#2) proportions without UNASSIGNED sites

proportions_assig<-matrix(0,dim(t_supmix)[[2]],5) ### activate this when calculating proportions including unknowns
for (i in 1:dim(t_supmix)[[2]]){
	important_col_ass<-important_col[,1:5]
	important_col_ass[i,]/apply(important_col_ass,1,sum)[i] -> proportions_assig[i,]
}
dim(proportions_assig)
head(proportions_assig)
summary(proportions_assig)

#Add columns of information Interval, chr, Cumulative,index

cbind(as.data.frame(chr[,c(1,4,6,7)]), as.data.frame(proportions_assig)) ->Pop_summary_assig
colnames(Pop_summary_assig)<-c("Interval","Chro","index","Cumulative","red","blue","green","orange","purple")
Pop_summary_assig ## remove yellow for just known things


# write the table with proportions of contribution of WBDC to eahc landrace population

output_contribution<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/Pop_summary/","Pop",p,"_w_unassigned_contributionFromWBDC.xls",sep="")
write.table(Pop_summary_unassi,output_contribution,quote=F,row.names=F,col.names=T,sep="\t")

output_contribution<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/Pop_summary/","Pop",p,"_all_assigned_contributionFromWBDC.xls",sep="")
write.table(Pop_summary_assig,output_contribution,quote=F,row.names=F,col.names=T,sep="\t")

#Get the average proportion of ancestry genomewide in each population

ancestry<-Pop_summary_unassi[,-c(1:4)]
average_ancestry<-apply(ancestry,2,mean)
SUMMARY_w_unassigned[p,]<-average_ancestry

ancestry2<-Pop_summary_assig[,-c(1:4)]
#Remove NaN rows
NaN_rows<-which(is.na(ancestry2[,1]))
ancestry2_noNA<-ancestry2[-c(NaN_rows),]
average_ancestry2<-apply(ancestry2_noNA,2,mean,is.na=F)
SUMMARY_all_assigned[p,]<-average_ancestry2

}

#Remove fourth row from Summary as we don't have a wild pop 4 (caspian sea) as a testing sample
SUMMARY_w_unassigned_plot<-SUMMARY_w_unassigned[-c(4),]
colnames(SUMMARY_w_unassigned_plot)<-c("N_Mesopotamia_training","S_Levant_training","Syrian_Desert_training","N_Levant_training","C_Asia_training","Unassigned")
row.names(SUMMARY_w_unassigned_plot)<-c("N_Mesopotamia","S_Levant","Syrian_Desert","N_Levant","C_Asia")
SUMMARY_w_unassigned_plot

write.table(SUMMARY_w_unassigned_plot,(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/Pop_summary/Genome_wide_prop_summary_w_unassigned.txt",sep="")),quote=F,row.names=T,col.names=T,sep="\t")

#for only assigned sites
SUMMARY_all_assigned_plot<-SUMMARY_all_assigned[-c(4),]
colnames(SUMMARY_all_assigned_plot)<-c("N_Mesopotamia_training","S_Levant_training","Syrian_Desert_training","N_Levant_training","C_Asia_training")
row.names(SUMMARY_all_assigned_plot)<-c("N_Mesopotamia","S_Levant","Syrian_Desert","N_Levant","C_Asia")
SUMMARY_all_assigned_plot

write.table(SUMMARY_all_assigned_plot,(paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Analysis/Pop_summary/Genome_wide_prop_summary_all_assigned.txt",sep="")),quote=F,row.names=T,col.names=T,sep="\t")

