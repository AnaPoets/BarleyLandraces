#!/bin/bash

#Create a directory for this rep

DIR=/home/SupportMix

#Create tpedData temp dir
mkdir $DIR/tpedData
#Change the window size you will be working. This to choose the right template configuration file
SNP=75SNP

#Set the directory to the directory containing your input files.
SNP_DIR=$DIR/75_SNPs

cd $SNP_DIR

#Change the coordinates from which to which run 
for subdir in {1..53}

do

cd $SNP_DIR/20_wild_testing_rep${subdir}

for i in *
do
gzip $i
done

#move gz files to tpedData to be pass to SupportMix

scp *.gz $DIR/tpedData

#Create a directory to keep original .gz files

mkdir tpedData
scp *.gz $SNP_DIR/20_wild_testing_rep${subdir}/tpedData

#Create directory to keep output 

mkdir Results
cd Results
mkdir chr1 chr2 chr3 chr4 chr5 chr6 chr7


#Go back to main directory

cd $DIR

for CHR in 1 2 3 4 5 6 7

do

SUPPORTMIX=/home/Scripts/SupportMixDistribution/Application/SupportMix

$SUPPORTMIX -c$CHR -C $DIR/Config_1896_crossval_chr${CHR}.cfg

mv outSupportMix_chr$CHR* $SNP_DIR/20_wild_testing_rep${subdir}/Results/chr$CHR

done

for OUTPUT in 1 2 3 4 5 6 7

do
scp $SNP_DIR/20_wild_testing_rep${subdir}/Results/chr$OUTPUT/outSupportMix_chr$OUTPUT.Probs.tped $SNP_DIR/20_wild_testing_rep${subdir}/Results/chr$OUTPUT/outSupportMix_chr$OUTPUT.tped $SNP_DIR/20_wild_testing_rep${subdir}

done

#Delete .gz files from the tpedData that holds the files just for the analysis

cd $DIR/tpedData
rm *gz

done
