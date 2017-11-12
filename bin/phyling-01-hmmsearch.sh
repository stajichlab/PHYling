#!/usr/bin/bash

LOG_FOLDER=logs
QUEUEING=parallel
CPU=4
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM variable"
 exit
fi

SCRIPT_DIR=$(dirname $0)
SUBJOB_SCRIPT=${SCRIPT_DIR}/phyling-01-hmmsearch-$QUEUEING.sh

JOBCOUNT=$(wc -l pepfile.lst)

if [[ $JOBCOUNT == 0 || $JOBCOUNT == "" ]]; then
    echo "No files to process, did you run PHYling init?"
    exit
fi

if [ $QUEUEING == "parallel" ]; then
    echo "run parallel job set"
    echo "parallel $SUBJOB_SCRIPT < pepfile.lst
elsif [ $QUEUEING == "slurm" ]; then
    echo "sbatch --array=1-$JOBCOUNT $SUBJOB_SCRIPT"
else
 echo "Just run in serial"
 for file in $(cat pepfile.lst)
 do
   echo "$SUBJOB_SCRIPT $file"
 done
fi

