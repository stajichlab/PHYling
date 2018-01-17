#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=combine
#SBATCH --time=2:00:00
#SBATCH --output=logs/combine.%A.log

if [ $MODULESHOME ]; then
    module switch python/3
fi

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM and other variables"
 exit
fi

if [ ! $PHYLING_DIR ]; then
    PHYLING_DIR=$(basename $0)../
fi

while getopts d:o:v:i:s:e:r:p: option
do
 case "${option}"
 in
 d) DIR=${OPTARG};;
 o) OUT=${OPTARG};;
 v) DEBUG=${OPTARG};;
 i) INCLUDE=${OPTARG};;
 s) SEED=${OPTARG};;
 e) EXT=$(OPTARG};;
 r) RAND=${OPTARG};;
 p) PARTITIONS=${OPTARG};;
 esac
done

ARGS=""
if [ $DEBUG ]; then
    ARGS="-v"
fi

if [ $INCLUDE ]; then
    ARGS+=" --include $INCLUDE"
fi

#MAYBE THSE SHOULD BE SET BY config.txt

if [ $SEED ]; then
    ARGS+=" --seed $SEED"
fi

if [ $RAND ]; then
    ARGS+=" --rand $RAND"
fi

if [ $EXT ]; then
    ARGS+=" --ext $EXT"
fi

if [ $PARTITIONS ]; then
    ARGS+=" -p $PARTITIONS"
fi

$PHYLING_DIR/util/combine_multiseq_aln.py -d $DIR -o $OUT $ARGS
