#Calculate the best value of K for Structure using Clusterdness

rm(list=ls())

str1<-read.table("~/Structure/output_k7/Landrace_K7_R1_q")
str2<-read.table("~/Structure/output_k7/Landrace_K7_R2_q")
str3<-read.table("~/Structure//output_k7/Landrace_Af_K7_R3_q")
str4<-read.table("~/Structure//output_k7/Landrace_Af_K7_R4_q")
str5<-read.table("~/Structure//output_k7/Landrace_Af_K7_R5_q")
str6<-read.table("~/Structure//output_k7/Landrace_Af_K7_R6_q")
str7<-read.table("~/Structure//output_k7/Landrace_Af_K7_R7_q")
str8<-read.table("~/Structure//output_k7/Landrace_Af_K7_R8_q")
str9<-read.table("~/Structure/output_k7/Landrace_Af_K7_R9_q")
str10<-read.table("~/Structure//output_k7/Landrace_Af_K7_R10_q")


#remove sample and fake column

str1<-str1[,-c(1,2)]
str2<-str2[,-c(1,2)]
str3<-str3[,-c(1,2)]
str4<-str4[,-c(1,2)]
str5<-str5[,-c(1,2)]
str6<-str6[,-c(1,2)]
str7<-str7[,-c(1,2)]
str8<-str8[,-c(1,2)]
str9<-str9[,-c(1,2)]
str10<-str10[,-c(1,2)]



Kvalue =7
Individuals = 803
#Determine G

#solve the square root

SumQ<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str1[i,k] -(1/Kvalue))^2) -> SumQ[i,k]
	apply(SumQ,1,sum)->SumQ2
	(SumQ2 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS
	
}
}

# Find G value

G_value4_1<-(1/Individuals)*sum(RESULTS)


#rep2

SumQ2<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS2<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str2[i,k] -(1/Kvalue))^2) -> SumQ2[i,k]
	apply(SumQ2,1,sum)->SumQ2_2
	(SumQ2_2 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS2
	
}
}

# Find G value

G_value4_2<-(1/Individuals)*sum(RESULTS2)

#Rep3
SumQ3<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS3<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str3[i,k] -(1/Kvalue))^2) -> SumQ3[i,k]
	apply(SumQ3,1,sum)->SumQ2_3
	(SumQ2_3 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS3
	
}
}

# Find G value

G_value4_3<-(1/Individuals)*sum(RESULTS3)

#Rep4
SumQ4<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS4<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str4[i,k] -(1/Kvalue))^2) -> SumQ4[i,k]
	apply(SumQ4,1,sum)->SumQ2_4
	(SumQ2_4 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS4
	
}
}

# Find G value

G_value4_4<-(1/Individuals)*sum(RESULTS4)

#rep5

SumQ5<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS5<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str5[i,k] -(1/Kvalue))^2) -> SumQ5[i,k]
	apply(SumQ5,1,sum)->SumQ2_5
	(SumQ2_5 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS5
	
}
}

# Find G value

G_value4_5<-(1/Individuals)*sum(RESULTS5)

#Rep6

SumQ6<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS6<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str6[i,k] -(1/Kvalue))^2) -> SumQ6[i,k]
	apply(SumQ6,1,sum)->SumQ2_6
	(SumQ2_6 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS6
	
}
}

# Find G value

G_value4_6<-(1/Individuals)*sum(RESULTS6)

#Rep7
SumQ7<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS7<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str7[i,k] -(1/Kvalue))^2) -> SumQ7[i,k]
	apply(SumQ7,1,sum)->SumQ2_7
	(SumQ2_7 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS7
	
}
}

# Find G value

G_value4_7<-(1/Individuals)*sum(RESULTS7)

#Rep8
SumQ8<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS8<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str8[i,k] -(1/Kvalue))^2) -> SumQ8[i,k]
	apply(SumQ8,1,sum)->SumQ2_8
	(SumQ2_8 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS8
	
}
}

# Find G value

G_value4_8<-(1/Individuals)*sum(RESULTS8)

#Re9

SumQ9<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS9<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str9[i,k] -(1/Kvalue))^2) -> SumQ9[i,k]
	apply(SumQ9,1,sum)->SumQ2_9
	(SumQ2_9 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS9
	
}
}

# Find G value

G_value4_9<-(1/Individuals)*sum(RESULTS9)

#Rep10

SumQ10<-matrix(nrow=Individuals,ncol=Kvalue)
RESULTS10<-NULL


for (i in 1:Individuals) {
	
	for (k in 1:Kvalue) {
	((str10[i,k] -(1/Kvalue))^2) -> SumQ10[i,k]
	apply(SumQ10,1,sum)->SumQ2_10
	(SumQ2_10 * (Kvalue/(Kvalue-1)))^0.5 ->RESULTS10
	
}
}

# Find G value

G_value4_10<-(1/Individuals)*sum(RESULTS10)


##Calculate the mean G value from 10 reps

mean(G_value4_1,G_value4_2,G_value4_3,G_value4_4,G_value4_5,G_value4_6,G_value4_7,G_value4_8,G_value4_9,G_value4_10) ->mean_G_value

cbind(G_value4_1,G_value4_2,G_value4_3,G_value4_4,G_value4_5,G_value4_6,G_value4_7,G_value4_8,G_value4_9,G_value4_10) ->Rep_value

##CALCULATE Confidence intervals
#import Clusterdness output
dat<-read.table("~/Clusteredness/10Rep_iSelect6152.txt",header=T)
head(dat)

STDEV<-dat$Stdev

CI<-NULL
for (i in 1:6) { qnorm(0.975)*(STDEV[i]/sqrt(10)) ->CI[i]}
CI

UpperCI<-dat$Average_10Runs + CI
LowerCI<-dat$Average_10Runs - CI

cbind(dat, as.data.frame(UpperCI),as.data.frame(LowerCI)) ->Results_plot
write.table(Results_plot,"~/Clusterdness/Results_plot.txt",row.names=F,col.names=T,sep="\t",quote=F)

