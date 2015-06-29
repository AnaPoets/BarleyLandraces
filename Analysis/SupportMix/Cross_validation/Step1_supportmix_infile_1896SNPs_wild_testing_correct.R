rm(list=ls())
#Define the NAME directory where output will be sent

#Create a directory for each window size analyzed
WindowsSize<-c("75_SNPs/")
DIRECTORY_windows<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,sep="")
dir.create(DIRECTORY_windows)

#Number of runs
for (D in 1:50) {

DirectoryName<-paste("20_wild_testing_rep",D,sep="")

#Create a folder for each rep to be done
DIRECTORY<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,sep="")

#Create the directory in the right place
dir.create(DIRECTORY)


#Create an directory called "Analysis"

ANALYSIS<-paste(DIRECTORY,"/Analysis",sep="")
dir.create(ANALYSIS)
	
SAMPLES<-read.table("~/Documents/SupporMix/SupportMix_wild_testing_1896/latlong_wild_double.txt",header=T)
head(SAMPLES)


#=============================
#Get .tped file for all the wild samples in the .tped file used for Supportmix with the landraces. Put them together
sup_lat_wild_pop1<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/input_HEADERS/ancestral_pop1.tped",header=T)
sup_lat_wild_pop2<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/input_HEADERS/ancestral_pop2.tped",header=T)
sup_lat_wild_pop3<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/input_HEADERS/ancestral_pop3.tped",header=T)
sup_lat_wild_pop5<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/input_HEADERS/ancestral_pop5.tped",header=T)
sup_lat_wild_pop6<-read.table("~/Documents/ASupporMix/Supportmix_barley_phased_1896/input_HEADERS/ancestral_pop6.tped",header=T)

#Separate chr_info
chr_info_land_wild<-sup_lat_wild_pop1[,c(1:4)]

geno_wild_pop1<-sup_lat_wild_pop1[,-c(1:4)]
geno_wild_pop2<-sup_lat_wild_pop2[,-c(1:4)]
geno_wild_pop3<-sup_lat_wild_pop3[,-c(1:4)]
geno_wild_pop5<-sup_lat_wild_pop5[,-c(1:4)]
geno_wild_pop6<-sup_lat_wild_pop6[,-c(1:4)]

geno_wild_ALL<-cbind(geno_wild_pop1, geno_wild_pop2, geno_wild_pop3, geno_wild_pop5, geno_wild_pop6)

#Get only wild samples, with names doubled

Wild_supmix<-geno_wild_ALL[,(colnames(geno_wild_ALL) %in% SAMPLES$Sample)]
dim(Wild_supmix)
#Wild genotypes and chr info

all_wild<-cbind(as.data.frame(chr_info_land_wild),as.data.frame(Wild_supmix))
#================================
#all_wild<-read.table("~/Dropbox/ANITA_301213/THESIS/Data/for_SupporMix/SupportMix_wild_testing_1896/all_ancestral_pop.tped.txt",header=T)

#Separate chr info

chr_info<-all_wild[,c(1:4)]

#All_wild doesn't have the caspian sea samples
#Select only a single copy of each sample name.Then double the names to select both haplotypes.
samples_single<-seq(1,nrow(SAMPLES),by=2)
SAMPLES_ODD<-SAMPLES[c(samples_single),]

pop1<-subset(SAMPLES_ODD, SAMPLES_ODD$pop6_hap == '1')
pop2<-subset(SAMPLES_ODD, SAMPLES_ODD$pop6_hap == '2')
pop3<-subset(SAMPLES_ODD, SAMPLES_ODD$pop6_hap == '3')
pop5<-subset(SAMPLES_ODD, SAMPLES_ODD$pop6_hap == '5')
pop6<-subset(SAMPLES_ODD, SAMPLES_ODD$pop6_hap == '6')

#Select subset of samples to be tested at each population
##POP1
#Select the columns with the wild samples to be tested, as a single representation per sample
sample_pop1<-sample(nrow(pop1),4,replace=F)
selected1<-pop1[c(sample_pop1),]
selected1_sample<-selected1$Sample

#Select the name of the columns for selected samples, to identify the double name as appear in all_wild
SAMPLED_pop1 <-data.frame()
for (i in 1:length(selected1_sample)){
	NAME<-selected1_sample[i]
SAMPLED_pop1<-as.matrix(append(SAMPLED_pop1,as.matrix(colnames(all_wild)[grep(paste("^",NAME,sep=""),colnames(all_wild))])))
}

##POP2

#Select the columns with the wild samples to be tested, as a single representation per sample
sample_pop2<-sample(nrow(pop2),4,replace=F)
selected2<-pop2[c(sample_pop2),]
selected2_sample<-selected2$Sample

#Select the name of the columns for selected samples, to identify the double name as appear in all_wild
SAMPLED_pop2 <-data.frame()
for (i in 1:length(selected2_sample)){
	NAME<-selected2_sample[i]
SAMPLED_pop2<-as.matrix(append(SAMPLED_pop2,as.matrix(colnames(all_wild)[grep(paste("^",NAME,sep=""),colnames(all_wild))])))
}

#POP3
#Select the columns with the wild samples to be tested, as a single representation per sample
sample_pop3<-sample(nrow(pop3),4,replace=F)
selected3<-pop3[c(sample_pop3),]
selected3_sample<-selected3$Sample

#Select the name of the columns for selected samples, to identify the double name as appear in all_wild
SAMPLED_pop3 <-data.frame()
for (i in 1:length(selected3_sample)){
	NAME<-selected3_sample[i]
SAMPLED_pop3<-as.matrix(append(SAMPLED_pop3,as.matrix(colnames(all_wild)[grep(paste("^",NAME,sep=""),colnames(all_wild))])))
}

#POP5

#Select the columns with the wild samples to be tested, as a single representation per sample
sample_pop5<-sample(nrow(pop5),4,replace=F)
selected5<-pop5[c(sample_pop5),]
selected5_sample<-selected5$Sample

#Select the name of the columns for selected samples, to identify the double name as appear in all_wild
SAMPLED_pop5 <-data.frame()
for (i in 1:length(selected5_sample)){
	NAME<-selected5_sample[i]
SAMPLED_pop5<-as.matrix(append(SAMPLED_pop5,as.matrix(colnames(all_wild)[grep(paste("^",NAME,sep=""),colnames(all_wild))])))
}

#POP6
#Select the columns with the wild samples to be tested, as a single representation per sample
sample_pop6<-sample(nrow(pop6),4,replace=F)
selected6<-pop6[c(sample_pop6),]
selected6_sample<-selected6$Sample

#Select the name of the columns for selected samples, to identify the double name as appear in all_wild
SAMPLED_pop6 <-data.frame()
for (i in 1:length(selected6_sample)){
	NAME<-selected6_sample[i]
SAMPLED_pop6<-as.matrix(append(SAMPLED_pop6,as.matrix(colnames(all_wild)[grep(paste("^",NAME,sep=""),colnames(all_wild))])))
}

#Put together all the samples selected. 20 samples = 40 haplotypes
Samples_Testing<-as.data.frame(rbind(SAMPLED_pop1, SAMPLED_pop2, SAMPLED_pop3, SAMPLED_pop5, SAMPLED_pop6))


##PREPARE TPED AND TFAM FILES FOR THE ADMIXED POPULATIONS
## Select all samples for admixed_wbdc_testing.tped file
#add c(1:4) this will call the columns for "Chro,SNP, EXTRA, gmap"
pop_admixed<-all_wild[,(colnames(all_wild) %in% Samples_Testing[,1])]
dim(pop_admixed)
#Add chr info
all_wild_admixed <-cbind(as.data.frame(chr_info),as.data.frame(pop_admixed))
all_wild_colnames<-colnames(all_wild_admixed)
all_wild_admixed[1:10,1:10]
write.table(all_wild_admixed,paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/admixed_wbdc_testing.tped",sep=""),quote=F,col.names=F,sep=" ",row.names=F)

#Print with headers
write.table(all_wild_admixed,paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/Header_admixed_wbdc_testing.tped",sep=""),quote=F,col.names=T,sep=" ",row.names=F)

#Construct tfam file for admixed (testing) individuals
#Remove chr info
admix_wild_colnames_only<-as.data.frame(all_wild_colnames[-c(1:4)])
#Select single names
SINGLE<-seq(1,nrow(admix_wild_colnames_only),2)
admix_wild_colnames_only_single<-as.data.frame(admix_wild_colnames_only[c(SINGLE),])

fake<-rep(0,length(admix_wild_colnames_only_single))

admixed_tfam_ready<-cbind(admix_wild_colnames_only_single, admix_wild_colnames_only_single,fake,fake,fake,fake)

write.table(admixed_tfam_ready,paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/admixed_wbdc_testing.tfam",sep=""),quote=F,col.names=F,sep="\t",row.names=F)


#PREPARE THE TPED AND TFAM FILE FOR EACH WILD (LEARNING) POPULATIONS
NAMES_ALL_WILD_LEARNING<-as.data.frame(setdiff(colnames(all_wild), colnames(all_wild_admixed)))

wild_learning<-all_wild[,(colnames(all_wild) %in% NAMES_ALL_WILD_LEARNING[,1])]

#CREATE TPED FILES FOR LEARNING SAMPLES BY POPULATION
#identify idenviduals in each population using names doubled

pop1_learning<-subset(SAMPLES,SAMPLES$pop6_hap == "1")
pop2_learning<-subset(SAMPLES,SAMPLES$pop6_hap == "2")
pop3_learning<-subset(SAMPLES,SAMPLES$pop6_hap == "3")
pop5_learning<-subset(SAMPLES,SAMPLES$pop6_hap == "5")
pop6_learning<-subset(SAMPLES,SAMPLES$pop6_hap == "6")

##POP1_LEARNING
pop_learning_1<-wild_learning[,(colnames(wild_learning) %in% pop1_learning$Sample)]
#add chromosome information
pop_learning_1_tped<-cbind(chr_info, pop_learning_1)

##POP2_LEARNING
pop_learning_2<-wild_learning[,(colnames(wild_learning) %in% pop2_learning$Sample)]
#add chromosome information
pop_learning_2_tped<-cbind(chr_info, pop_learning_2)

#POP3_LEARNING
pop_learning_3<-wild_learning[,(colnames(wild_learning) %in% pop3_learning$Sample)]
#add chromosome information
pop_learning_3_tped<-cbind(chr_info, pop_learning_3)

#POP5_LEARNING
pop_learning_5<-wild_learning[,(colnames(wild_learning) %in% pop5_learning$Sample)]
#add chromosome information
pop_learning_5_tped<-cbind(chr_info, pop_learning_5)

#POP6_LEARNING
pop_learning_6<-wild_learning[,(colnames(wild_learning) %in% pop6_learning$Sample)]
#add chromosome information
pop_learning_6_tped<-cbind(chr_info, pop_learning_6)

#WRITE TPED ANCESTRAL FILES

for (i in c(1,2,3,5,6)) {outfile <-paste("ancestral_pop",i,".tped",sep="")
	directory<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/",outfile,sep="")
	write.table(get(paste("pop_learning_",i,"_tped",sep="")),directory,quote=F,col.names=F,sep=" ",row.names=F)}

##PREPARE TFAM FILES FOR ANCESTRAL POPS
#Construct tfam file for admixed (testing) individuals
#Remove chr info

#tfam wild1
wild1_learning_colnames<-colnames(pop_learning_1_tped)

wild1_learning_names<-as.data.frame(wild1_learning_colnames[-c(1:4)])
SINGLE_w<-seq(1,nrow(wild1_learning_names),2)
wild1_colnames_only_single<-as.data.frame(wild1_learning_names[c(SINGLE_w),])

fake<-rep(0,length(wild1_colnames_only_single))

wild1_tfam_ready<-cbind(wild1_colnames_only_single, wild1_colnames_only_single,fake,fake,fake,fake)

#tfam wild2 ancestral
wild2_learning_colnames<-colnames(pop_learning_2_tped)

wild2_learning_names<-as.data.frame(wild2_learning_colnames[-c(1:4)])
#Select single names
SINGLE_w<-seq(1,nrow(wild2_learning_names),2)
wild2_colnames_only_single<-as.data.frame(wild2_learning_names[c(SINGLE_w),])

fake<-rep(0,length(wild2_colnames_only_single))

wild2_tfam_ready<-cbind(wild2_colnames_only_single, wild2_colnames_only_single,fake,fake,fake,fake)


#tfam wild3 ancestral
wild3_learning_colnames<-colnames(pop_learning_3_tped)

wild3_learning_names<-as.data.frame(wild3_learning_colnames[-c(1:4)])
#Select single names
SINGLE_w<-seq(1,nrow(wild3_learning_names),2)
wild3_colnames_only_single<-as.data.frame(wild3_learning_names[c(SINGLE_w),])

fake<-rep(0,length(wild3_colnames_only_single))

wild3_tfam_ready<-cbind(wild3_colnames_only_single, wild3_colnames_only_single,fake,fake,fake,fake)

#tfam wild25 ancestral
wild5_learning_colnames<-colnames(pop_learning_5_tped)

wild5_learning_names<-as.data.frame(wild5_learning_colnames[-c(1:4)])
#Select single names
SINGLE_w<-seq(1,nrow(wild5_learning_names),2)
wild5_colnames_only_single<-as.data.frame(wild5_learning_names[c(SINGLE_w),])

fake<-rep(0,length(wild5_colnames_only_single))

wild5_tfam_ready<-cbind(wild5_colnames_only_single, wild5_colnames_only_single,fake,fake,fake,fake)

#tfam wild6 ancestral
wild6_learning_colnames<-colnames(pop_learning_6_tped)

wild6_learning_names<-as.data.frame(wild6_learning_colnames[-c(1:4)])
#Select single names
SINGLE_w<-seq(1,nrow(wild6_learning_names),2)
wild6_colnames_only_single<-as.data.frame(wild6_learning_names[c(SINGLE_w),])

fake<-rep(0,length(wild6_colnames_only_single))

wild6_tfam_ready<-cbind(wild6_colnames_only_single, wild6_colnames_only_single,fake,fake,fake,fake)


#write.table(pop6_tfam_ready,"~/Dropbox/ANITA_301213/THESIS/Data/for_SupporMix/SupportMix_wild_testing_1896/20_wild_testing_rep7/ancestral_pop6.tfam",quote=F,col.names=F,sep="\t",row.names=F)
for (i in c(1,2,3,5,6)) {outfile <-paste("ancestral_pop",i,".tfam",sep="")
	directory<-paste("~/Documents/SupporMix/SupportMix_wild_testing_1896/",WindowsSize,DirectoryName,"/",outfile,sep="")
	write.table(get(paste("wild",i,"_tfam_ready",sep="")),directory,quote=F,col.names=F,sep="\t",row.names=F)}

}

#Run All_WildTesting_1896.sh. The contig files for all the chromosomes have to be in the same 75_SNPs directory.