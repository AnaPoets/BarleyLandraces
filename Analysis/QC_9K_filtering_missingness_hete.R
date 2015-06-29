# remove all existing files
rm(list = ls())

#input file rows are Marker names and colums are samples names
alch<-read.table('/Users/gonzales/Documents/Anita/candidate_genes_barley/alch/Alch9k_matrix.xls"',header=T)

#have a list of SNPs names (No "SNP" header needed)
alch[,1] ->snp
as.data.frame(snp)->snp
# remove row of marker names otherwise they will appear as one more argument instead of a header
alch <- alch[,-1]

#convert input file (data frame) to a matrix
as.matrix(alch)->alch

#add the snp's names as row.names
row.names(alch)<-snp[,1]
#turn table 
t(alch)->talch
#convert file again to data.frame
as.data.frame(talch)->alch

#Order the file based on sample names (first column)
alch <- alch[order(row.names(alch)),]


# size of data set

dim(alch)

#alch <- rownames(alch[,1])

# routine to identify SNPs with missing data, shown as '--'
missing <- function(dat) { 
dat <- na.omit(dat)
miss <- dat[dat == "--"]
miss <- length(miss)
dat_size <- length(dat)
miss_locus <- miss/dat_size
return(miss_locus)
}

# routine to identify monomorphic SNPs
# actually identifies minor allele frequency (MAF)
mono <- function(dat) {
#dat_size <- length(na.omit(dat))
AA <- dat[dat == "AA"]
AA <- length(AA)
BB <- dat[dat == "BB"]
BB <- length(BB)
smaller <- min(c(AA,BB))
return(smaller)
}

# routine to identify heterozygous SNPs
hets <- function(dat) {
dat_size <- length(na.omit(dat))
het <- dat[dat == "AB"]
het <- length(het)
hets_locus <- het/dat_size
return(hets_locus)
}

# all of the functions below are run using apply
# remove SNPs with missing data ≥ 10% 
# remove SNPs with observed heterozygosity > 10%
# remove SNPs with diversity = 0%; MAF = 0
miss_rm <- as.vector(which(apply(alch,2,missing) >= 0.10))
het_rm <- as.vector(which(apply(alch,2,hets) >= 0.10))
mono_rm <- as.vector(which(apply(alch,2,mono) == 0))
#miss_rm <- as.vector(which(apply(alch,2,missing) >= 0.25))
#het_rm <- as.vector(which(apply(alch,2,hets) >= 0.09))
#mono_rm <- as.vector(which(apply(alch,2,mono) == 0))

alch <- alch[,c(-miss_rm,-het_rm,-mono_rm)]

# remove lines with missing data ≥ 10% 
# remove lines with observed heterozygosity > 10%
miss_rm_accessions <- as.vector(which(apply(alch,1,missing) >= 0.10))
het_rm_accessions <- as.vector(which(apply(alch,1,hets) >= 0.10))

# now write NA values where there 
# het_count <- apply(alch,1,hets)

alch <- alch[c(-miss_rm_accessions,-het_rm_accessions),]

dim(alch)

#turn the table back to original format col.names=samples, row.names=markers

t(alch)->alch

#Write new table into a file
write.table(alch, file="~/Desktop/alch_95_matrix_good_Filtered.xls",quote=F,sep="\t",row.names=T,col.names=T)
