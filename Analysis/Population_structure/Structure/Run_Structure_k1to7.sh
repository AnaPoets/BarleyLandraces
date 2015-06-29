#!/bin/bash -l
# written by Peter L. Morrell, April 2011, St. Paul, MN 
#PBS -l walltime=70:00:00,mem=4gb,nodes=4:ppn=1
#PBS -m abe -M agonzale@umn.edu

# the full path to user space can vary
# the line below makes is easy to set
Working_dir=/home/morrellp/gonzales/Scripts
# directory containing Structure application
DIR=$Working_dir/Structure2.3.3/iSelect_6152_all_SNPs
# full path to application
PROGRAM=/home/morrellp/gonzales/Scripts/Structure2.3.3/console/structure
# Structure data file
INPUT1=$DIR/Str_input_1652.txt
# Structure main parameters file
INPUT2=$DIR/Landr_mainparams
# Structure extra parameters file
INPUT3=$DIR/Landr_extraparams
# output file; name will be appended 
OUTPUT=$DIR/output/Landrace_results_

REPS=10
KMIN=1
KMAX=7

for ((j=$KMIN;j<=$KMAX;j++))
do

for ((i=1;i<=$REPS;i++))
do
    cd $DIR
    $PROGRAM -i $INPUT1 -m $INPUT2 -e $INPUT3 -K $j -o ${OUTPUT}_K${j}_R${i}
done
done
