#!/usr/bin/bash

LOG_FOLDER=logs
QUEUEING=parallel
JOBCPU=1
TOTALCPU=4
PEPEXT=aa.fasta
ALN_OUTDIR=aln
ALNTOOL=hmmalign
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM and other variables"
 exit
fi

JOBCPU=1 # force this to 1 since all of these steps are single threaded

if [ ! $HMM ]; then
 echo "need to a config file to set the HMM folder name"
fi

if [ $1 ]; then
    ALNTOOL=$1
fi
SCRIPT_DIR=$(dirname $0)
SUBJOB_SCRIPT=${SCRIPT_DIR}/phyling-03-aln-$ALNTOOL.sh

if [ $QUEUEING == "parallel" ]; then
    JOBPARALLEL=$(expr $TOTALCPU / $JOBCPU)
    echo "Run parallel job $ALNTOOL"
    #echo "$JOBPARALLEL $SUBJOB_SCRIPT"
    echo "ls $ALN_OUTDIR/$HMM/*.$PEPEXT | parallel -j $JOBPARALLEL $SUBJOB_SCRIPT"
elif [ $QUEUEING == "slurm" ]; then
    echo "Run srun job hmmsearch"
    for INPUTFILE in $(ls $ALN_OUTDIR/$HMM/*.$PEPEXT)
    do
	sbatch --ntasks $JOBCPU --nodes 1 --export=IN=$INPUTFILE $SUBJOB_SCRIPT
    done
else
 echo "Run in serial"
 for file in $(ls $ALN_OUTDIR/$HMM/*.$PEPEXT)
 do
   $SUBJOB_SCRIPT $file
 done
fi
