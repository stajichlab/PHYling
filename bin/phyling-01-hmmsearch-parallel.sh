#!/usr/bin/env bash

#SBATCH --nodes 1 --ntasks 4 --mem 2G --out=logs/hmmsearch.%A.log

if [ $MODULESHOME ]; then
    module load hmmer/3
fi

PEPDIR=pep
PEPEXT=aa.fasta
HMM_FOLDER=HMM
HMMSEARCH_CUTOFF=1e-10
HMMSEARCH_OUT=search
LOG_FOLDER=logs
CPU=2

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM variable"
 exit
fi

if [ ! $PHYLING_DIR ]; then
    echo "need a PHYLING_DIR defined in config.txt or env variable"
    exit
fi

if [ ! $HMM ]; then
 echo "need to a config file to set the HMM folder name"
fi

if [ ! $IN ]; then
    IN=$1
    if [ ! $IN ]; then
	echo "NO INPUT GENOME FILE"
	exit
    fi
fi

IN=$(basename $IN)
NM=$(basename $IN .$PEPEXT)
echo "g=$IN NM=$NM"

MARKERS=${HMM_FOLDER}/$HMM/markers_3.hmmb
OUT=$HMMSEARCH_OUT/$HMM
mkdir -p $OUT

# number of processors to use set by 
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
# convention is they all end in .aa.fasta - change this if not or make a variable

if [[ ! -f "$OUT/$NM.domtbl" ||
      $PEPDIR/$IN -nt $OUT/$NM.domtbl ]]; then
    hmmsearch --cpu $CPU -E $HMMSEARCH_CUTOFF --domtblout $OUT/$NM.domtbl $MARKERS $PEPDIR/$IN >& $OUT/$NM.log
else
 echo "skipping $NM - has already run"
fi

if [ ! -f $OUT/$NM.best ]; then
    ${PHYLING_DIR}/util/get_best_hmmtbl.py -c $HMMSEARCH_CUTOFF --input $OUT/$NM.domtbl > $OUT/$NM.best
fi
