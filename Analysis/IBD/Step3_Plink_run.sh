#!/bin/bash

#PBS -l mem=1000mb,nodes=1:ppn=1,walltime=72:00:00
#PBS -m abe
#PBS -M username@email.edu

module load plink
DIR=/home/PLINK/Sep_30SNP
cat List_gral.txt | while read line
do

SEGMENT=`echo $line`

plink --ped $DIR/${SEGMENT}.ped --map $DIR/${SEGMENT}.map --genome gz --out $DIR/output/${SEGMENT}out

done


