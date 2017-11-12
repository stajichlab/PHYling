#!/usr/bin/env bash

if [ $MODULESHOME != "" ]; then
    module load hmmer/3
fi
PEPDIR=pep
HMM_FOLDER=HMM
HMMSEARCH_CUTOFF=1e-10
HMMSEARCH_OUT=search
LOG_FOLDER=logs
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM variable"
 exit
fi

if [ ! $HMM ]; then
 echo "need to a config file to set the HMM folder name"
fi

G=$1
if [ ! $G ]; then
    echo "NO INPUT GENOME FILE"
    exit
fi

NM=`basename $G .aa.fasta`
echo "g=$G"

MARKERS=${HMM_FOLDER}/$HMM/markers_3.hmmb
OUT=$HMMSEARCH_OUT/$HMM
mkdir -p $OUT

# number of processors to use set by 
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
# convention is they all end in .aa.fasta - change this if not or make a variable

if [[ ! -f "$OUT/$NM.domtbl" || $PEPDIR/$G -nt $OUT/$NM.domtbl ]]; then
 echo "hmmsearch --cpu $CPU -E $CUTOFF --domtblout $OUT/$NM.domtbl $MARKERS $PEPDIR/$G >& $OUT/$NM.log"
else
 echo "skipping $NM - has already run"
fi
