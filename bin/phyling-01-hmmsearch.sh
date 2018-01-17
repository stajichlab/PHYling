#!/usr/bin/env bash

QUERYDBS=pepfile.lst
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
SUBJOB_SCRIPT=${SCRIPT_DIR}/phyling-01-hmmsearch-run.sh

if [ $QUEUEING == "parallel" ]; then
    JOBPARALLEL=$(expr $TOTALCPU / $JOBCPU)
    echo "Run parallel job hmmsearch"
    #echo "$JOBPARALLEL $SUBJOB_SCRIPT"
    parallel -j $JOBPARALLEL $SUBJOB_SCRIPT < $QUERYDBS
elif [ $QUEUEING == "slurm" ]; then
    QUEUECMD=""
    if [ $QUEUE ]; then
	QUEUECMD="-p $QUEUE"
    fi
    FILECOUNT=$(wc -l $QUERYDBS | awk '{print $1}')
    sbatch --ntasks $JOBCPU --nodes 1 $QUEUECMD --export=PHYLING_DIR=$0 \
	--array=1-${FILECOUNT} $SUBJOB_SCRIPT
else
 echo "Run in serial"
 for file in $(cat $QUERYDBS)
 do
     $SUBJOB_SCRIPT $file
 done
fi
