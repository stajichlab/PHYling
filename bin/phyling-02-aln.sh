#!/usr/bin/env bash

#defaults which are overridden config.txt if set
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
 exit
fi

if [ ! $ALNFILES ]; then
    ALNFILES=alnlist.$HMM
fi

if [[ ! -f $ALNFILES || $CLEAN ]]; then
    find $ALN_OUTDIR/$HMM -name "*.$PEPEXT" -exec basename {} \; > $ALNFILES
fi

if [ $1 ]; then
    ALNTOOL=$1
fi

SCRIPT_DIR=$(dirname $0)
SUBJOB_SCRIPT=${SCRIPT_DIR}/phyling-03-aln-$ALNTOOL.sh

if [ $QUEUEING == "parallel" ]; then
    JOBPARALLEL=$(expr $TOTALCPU "/" $JOBCPU)
    echo "Run parallel job $ALNTOOL"
    #echo "$JOBPARALLEL $SUBJOB_SCRIPT"
    parallel -j $JOBPARALLEL $SUBJOB_SCRIPT < $ALNFILES
    # do combine
    
elif [ $QUEUEING == "slurm" ]; then
    echo "Run srun job hmmsearch"
    #parallel -j $JOBPARALLEL srun --ntasks 1 --nodes 1 --export=IN={} $SUBJOB_SCRIPT 
    QUEUECMD=""
    if [ $QUEUE ]; then
	QUEUECMD="-p $QUEUE"
    fi
    ALNCT=$(wc -l $ALNFILES | awk '{print $1}')
    PHYLING_DIR=$(dirname $0)
    echo "PHYLING_DIR is $PHYLING_DIR"
    submid=$(sbatch --ntasks $JOBCPU --nodes 1 $QUEUECMD --export=PHYLING_DIR=$PHYLING_DIR \
	--array=1-${ALNCT} $SUBJOB_SCRIPT)

else
 echo "Run in serial"
 for file in $( find $ALN_OUTDIR/$HMM -name "*.$PEPEXT")
 do
   $SUBJOB_SCRIPT $file
 done
fi
