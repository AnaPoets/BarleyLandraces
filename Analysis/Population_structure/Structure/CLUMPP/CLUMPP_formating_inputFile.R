#Using 'cat' put together all the results for each K value.
# Then upload the resulting file "all_k_land_af_eu_as_reps.txt"

rm(list=ls())

for (i in 1:7){
all_k<-read.table(paste("~/Documents/land_iSelect_6135_k1-k7_rep/output_k",i,"/Landrace_Af_eu_as_results_iSelect6152_K4__all.txt",sep=""))
head(all_k)
dim(all_k)->dimentions_k
#Create a column of (1)
extra<-rep("(1)",dimentions_k[1])
#Create a column of pop# (all same number)
pop<-rep("8",dimentions_k[1])
semicolon<-rep(":",dimentions_k[1])
#devide in 10 groups for the 10 reps
number_sample<-rep(1:(dimentions_k[1]/10),10)

#Reformat table for CLUMPP
cbind(number_sample,number_sample,extra,pop,semicolon,all_k[,3:dimentions_k[[2]]]) ->new_all_k
head(new_all_k)

write.table(new_all_k,paste("~/Documents/Structure_output_reps/land_iSelect_6152_k1-k7_rep/output_k",i,"/CLUMPP_input_k4_all.txt",sep=""),quote=F, row.names=F,col.names=F,sep="\t")
}


