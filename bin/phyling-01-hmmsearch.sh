#!/usr/bin/env bash

LISTFILE=pepfile.lst
LOG_FOLDER=logs
QUEUEING=parallel
JOBCPU=2
TOTALCPU=4
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM variable"
 exit
fi

SCRIPT_DIR=$(dirname $0)
SUBJOB_SCRIPT=${SCRIPT_DIR}/phyling-01-hmmsearch-$QUEUEING.sh

if [ $QUEUEING == "parallel" ]; then
    JOBPARALLEL=$(expr $TOTALCPU / $JOBCPU)
    echo "Run parallel job hmmsearch"
    #echo "$JOBPARALLEL $SUBJOB_SCRIPT"
    parallel -j $JOBPARALLEL $SUBJOB_SCRIPT < $LISTFILE
elif [ $QUEUEING == "slurm" ]; then
    echo "Run srun job hmmsearch"
    for INPUTFILE in $(cat $LISTFILE)
    do
	sbatch --ntasks $JOBCPU --nodes 1 --export=PHYLING_DIR=$0 --export=IN=$INPUTFILE $SUBJOB_SCRIPT
    done
else
 echo "Run in serial"
 for file in $(cat $LISTFILE)
 do
   $SUBJOB_SCRIPT $file
 done
fi
