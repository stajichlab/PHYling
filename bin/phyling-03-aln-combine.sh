#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=combine
#SBATCH --time=2:00:00
#SBATCH --output=logs/combine.%A.log

if [ ! -z "$MODULESHOME" ]; then
    module unload miniconda2
    module load miniconda3 # UCR python3 enforcement
fi

if [ -e config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM and other variables"
 exit
fi

if [ -z "$PHYLING_DIR" ]; then
    PHYLING_DIR=$(basename $0)../
fi
DIR=$ALN_OUTDIR/$HMM


while getopts d:o:v:i:s:e:r:t:p:x: option
do
 case "${option}"
 in
 d) DIR=${OPTARG} ;;
 o) OUT=${OPTARG} ;;
 v) DEBUG=${OPTARG} ;;
 i) INCLUDE=${OPTARG} ;;
 s) SEED=${OPTARG} ;;
 e) EXT=${OPTARG} ;;
 r) RAND=${OPTARG} ;;
 t) TYPE=${OPTARG} ;;
 p) PARTITIONS=${OPTARG} ;;
 x) EXPECTED=${OPTARG} ;;
 esac
done

F=$EXPECTED
if [ -z "$F" ]; then
    F=$LISTFILE
fi
echo "reading $F to get count"
COUNT=$(wc -l $F | awk '{print $1}')
if [[ ! $COUNT || $COUNT == "0" ]]; then
    echo "NO QUERYDBS variable on config.txt?"
    COUNT="XX"
fi
OUTFILE=$PREFIX.${COUNT}_taxa.${HMM}.${TYPE}.fasaln
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

if [ $EXPECTED ]; then
 ARGS+=" --expected $EXPECTED"
fi
echo "$PHYLING_DIR/util/combine_multiseq_aln.py -d $DIR -o $OUTFILE $ARGS"
$PHYLING_DIR/util/combine_multiseq_aln.py -d $DIR -o $OUTFILE $ARGS
