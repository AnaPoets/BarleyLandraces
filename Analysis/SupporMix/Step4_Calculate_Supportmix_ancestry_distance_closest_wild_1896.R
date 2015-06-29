rm(list=ls())
support<-read.table("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/gmap_and_suport_complete.xls",header=T)

support[1:10,1:15]
#Remove SNP info

support<-support[,-c(1:12 )]

support[1:10,1:10] 
t(support)->t_support
t_support[1:10,]
table(t_support[2,])
dim(t_support)

#remove extra
t_support<-t_support[-1,]
#t_support<-t_support[1:10,]

matrix(0,dim(t_support)[[1]],7)->my_results

for (i in 1:dim(t_support)[1]){
	table(t_support[i,]) ->table_R
	as.data.frame(table_R)->my_table
	intersect(c(1,2,3,5,6,7),my_table[,1]) ->table_data
	as.data.frame(table_data)->items_value
	
	for (j in 1:dim(items_value)[1]) {
	items_value[j,1] ->position
	as.numeric(as.character(position))->position2
	     my_table[j,2] ->my_results[i,position2]}
	}
	
	my_results[1:10,]
	dim(my_results)
	colnames(my_results)<-c("1","2","3","4","5","6","7")
	important_col<-my_results[,c(1,2,3,5,6)]
	important_col[1:10,]
	
#Add sample names

row.names(important_col)<-row.names(t_support)

important_col[1:10,]
colnames(important_col)<-c("Northern_Mesopotamia","Southern_Levant","Syrian_Desert","Northern_levant","Central_Asia")
important_col<-as.data.frame(important_col)
#Calculate average proportion at each sample (23=regions with no missing data)

proportion_ancestry_samples <-important_col/sum(important_col[1,])
proportion_ancestry_samples[1:10,]

if (length(which(is.na(proportion_ancestry_samples))) ==0) {

#get only one sample (allele)

UNIQUE<-seq(1,nrow(proportion_ancestry_samples),2)
proportion_ancestry<-proportion_ancestry_samples[UNIQUE,]
dim(proportion_ancestry)
#Add geographic location

latlong<-read.table("~/Documents/github/BarleyLandraces/Datasets/Samples_latlong.txt",header=T)
head(latlong)
latlong<-latlong[,-c(2,3,4,5,6,7,8)]

dim(latlong)

#Select only landraces
latlong_land<-latlong[(latlong[,1]%in% row.names(proportion_ancestry)),]
#sort by sample name

latlong_l_or<-latlong_land[order(latlong_land[,1]),]

#Attach latlong to proportion ancestry

prop_latlong<-cbind(as.data.frame(latlong_l_or),as.data.frame(proportion_ancestry))

######
# Calculate the geodesic distance between two points specified by degrees (DD) latitude/longitude using
# Haversine formula (hf), Spherical Law of Cosines (slc) and Vincenty inverse formula for ellipsoids (vif)
library(pracma)
# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

gcd <- function(long1, lat1, long2, lat2) {
 
  # Convert degrees to radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
 
  return(haversine = gcd.hf(long1, lat1, long2, lat2)                 
	
)
}


#####

#Calculating distance away from Central Asia

wild<-read.table("~/Documents/github/BarleyLandraces/Datasets/WBDC/WBDC_assig_k6_hap.txt",header=T)
wild[1:10,]

#Separate wild Populations
w_n_mesopotamia<-subset(wild,wild$pop6_hap == '1')
w_s_levant<-subset(wild,wild$pop6_hap == '2')
w_s_desert<-subset(wild,wild$pop6_hap == '3')
w_n_levant<-subset(wild,wild$pop6_hap == '5')
w_asia<-subset(wild,wild$pop6_hap == '6')

#####Choose the wild population to use as reference
POPULATIONS<-c("w_n_mesopotamia","w_s_levant","w_s_desert","w_n_levant", "w_asia")
Pop_names<-c("N_mesopotamia","S_levant","S_desert","N_levant","C_asia")
for (p in 1:length(POPULATIONS)){
population<-get(POPULATIONS[p])

Distance<-NULL
land_distance<-NULL

#Compare one landrace to each wild to determine shortest distance
for (i in 1:(dim(prop_latlong)[1])) {
	for (w in 1:(dim(population)[1])) {
	
	Distance[w]<-gcd(prop_latlong[i,3],prop_latlong[i,2],population[w,4],population[w,3])
	land_distance[i]<-min(Distance)
	}
	}
#} else print("Error! there are samples with NaN == no assignment")
prop_latlong$Distance<-land_distance
length(land_distance)
head(prop_latlong)
dim(prop_latlong)

write.table(prop_latlong,paste("~/Documents/SupporMix/Supportmix_barley_phased_1896/Analysis/75_windows/Ancestry_vs_distance/From_", Pop_names[p],".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
}
