#First run Supportmix_out_output_1896_step2.R

# The following code will calculate the proportion of each Wild population for each landrace population by each chromosomal segment. And PLOT
rm(list=ls())

#Set working directory's name

supmix<-read.table("~/DocumentsSupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/gmap_and_suport_complete.xls",header=T)
dim(supmix)
supmix[1:10,1:14]
#Select the chrom windows information (cM, snp name, etc)
chr<-supmix[,1:13]
supmix_all_samples<-supmix[,-c(1:13)]


#call the population assignement according to Structure k=4.
#Separate wild barley by populations
str_land<-read.table("~/Documents/land_assignment_k4_iSelect6135.txt",header=T)
head(str_land)

#make another str_wild file adding the '_2' to each sample name. Then combine the files to later separate all the wild samples.
str_land_2<-cbind(paste(str_land[,1],'_2',sep=""),as.data.frame(str_land[,-c(1)]))
colnames(str_land_2)<-c("Samples" ,"LATITUDE", "LONGITUDE" ,"Populations")

assi_str_all<-rbind(as.data.frame(str_land),as.data.frame(str_land_2))
head(assi_str_all)


#Create a subdirectory to keep the output from proportions
dir.create(paste("~/Dropbox/ANITA_301213/THESIS/Data/for_SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary",sep=""))

# Select samples names for each population

pop1<-subset(assi_str_all, assi_str_all $Populations == "1")
pop2<-subset(assi_str_all, assi_str_all $Populations == "2")
pop3<-subset(assi_str_all, assi_str_all $Populations == "3")
pop4<-subset(assi_str_all, assi_str_all $Populations == "4")
# Separate resutls from supportMix for each population

supmix[,(colnames(supmix) %in% pop1$Samples)]->supmix1
supmix[,(colnames(supmix) %in% pop2$Samples)]->supmix2
supmix[,(colnames(supmix) %in% pop3$Samples)]->supmix3
supmix[,(colnames(supmix) %in% pop4$Samples)]->supmix4


#Paste the chr information to each supmix population

cbind(as.data.frame(chr),as.data.frame(supmix1)) -> chr_supmix1
cbind(as.data.frame(chr),as.data.frame(supmix2)) -> chr_supmix2
cbind(as.data.frame(chr),as.data.frame(supmix3)) -> chr_supmix3
cbind(as.data.frame(chr),as.data.frame(supmix4)) -> chr_supmix4

##Calculate contribution from wild to each landrace population
#for each population (p)analyze the proportions
for (p in 1:4){

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
	#WithOUT unassigned populations
	important_col<-my_results[,c(1,2,3,5,6)]
	

#Calculate in proportions what is the contribution of each Wild (7 columns) to the landraces population at each segment

matrix(0,dim(t_supmix)[[2]],5)->proportions ### activate this when calculating proportions including unknowns
for (i in 1:dim(t_supmix)[[2]]){
	important_col[i,]/apply(important_col,1,sum)[i] ->proportions[i,]
}

#Remove rows of NaN due to having only Unassigned data
NaN_remove<-which(is.na(proportions[,1]))
proportions<-proportions[-(NaN_remove),]

chr_removed<-chr[-(NaN_remove),]
dim(proportions)
head(proportions)
summary(proportions)
#Add columns of information Interval, chr, Cumulative,index

cbind(as.data.frame(chr_removed[,c(1,4,6,7)]), as.data.frame(proportions)) ->Pop_summary
colnames(Pop_summary)<-c("Interval","Chro","index","Cumulative","red","blue","green","orange","purple")
Pop_summary


output_contribution<-paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Pop_summary/","Pop",p,"_contributionFromWBDC.xls",sep="")
write.table(Pop_summary,output_contribution,quote=F,row.names=F,col.names=T,sep="\t")
}

#Continue with  Step4_Plot_SupportMix_deltaAncestry.R