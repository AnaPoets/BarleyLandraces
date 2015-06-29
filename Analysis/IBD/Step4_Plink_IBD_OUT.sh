#!/bin/bash

#This script will take the output files from PLINK and select those shared haplotypes with high similarity "1.0000". It will retain only
# those pairs that involve a sample from the wild with a sample from the landraces.
# The output will be sent to a new directory "Analysis"
DIR=~/Documents/PLINK/Sep_30SNP/output

cd $DIR
mkdir Analysis

cat $DIR/List_gral.txt | while read line
do

SEGMENT=`echo $line`

gzcat ${SEGMENT}out.genome.gz | awk '$10 == "1.0000" {print $0}'|awk '{if ($2 ~ /^WBDC/ && $4 !~ /^WBDC/) {print $0} else if ($2 !~ /^WBDC/ && $4 ~ /^WBDC/) {print $0}} ' >$DIR/Analysis/${SEGMENT}_perfectMatch.txt
done
