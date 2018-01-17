#!/usr/bin/env bash

#defaults which are overridden config.txt if set
LOG_FOLDER=logs
QUEUEING=parallel
JOBCPU=1
TOTALCPU=4
OUTPEPEXT=aa.fa
OUTCDSEXT=cds.fa
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

while getopts t:c:f: option
do
 case "${option}"
 in
 t) ALNTOOL=${OPTARG};;
 c) CLEAN=${OPTARG};;
 f) FORCE=${OPTARG};;
 esac
done

if [[ ! -f $ALNFILES || $CLEAN || $FORCE ]]; then
    find $ALN_OUTDIR/$HMM -name "*.$OUTPEPEXT" > $ALNFILES
fi

SCRIPT_DIR=$(dirname $0)
SUBJOB_SCRIPT=${SCRIPT_DIR}/phyling-02-aln-$ALNTOOL.sh

if [ $QUEUEING == "parallel" ]; then
    JOBPARALLEL=$(expr $TOTALCPU "/" $JOBCPU)
    echo "Run parallel job $ALNTOOL"
    #echo "$JOBPARALLEL $SUBJOB_SCRIPT"
    parallel -j $JOBPARALLEL $SUBJOB_SCRIPT -f $FORCE -i {} < $ALNFILES
    # do combine here
    
elif [ $QUEUEING == "slurm" ]; then
    QUEUECMD=""
    if [ $QUEUE ]; then
	QUEUECMD="-p $QUEUE"
    fi
    ALNCT=$(wc -l $ALNFILES | awk '{print $1}')
    PHYLING_DIR=$(dirname $0)
    echo "PHYLING_DIR is $PHYLING_DIR"
    submitid=$(sbatch --ntasks $JOBCPU --nodes 1 $QUEUECMD \
	--export=PHYLING_DIR=$PHYLING_DIR --export=FORCE=$FORCE \
	--array=1-${ALNCT} $SUBJOB_SCRIPT)
    
    # do combining as another slurm job
else
 echo "Run in serial"
 for file in $( find $ALN_OUTDIR/$HMM -name "*.$OUTPEPEXT")
 do
   $SUBJOB_SCRIPT -f $FORCE $file
 done
 
 # do combine here (how to avoid duplicate code with the parallel stuff above

fi

