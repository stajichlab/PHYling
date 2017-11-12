#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=hmmsearch
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmsearch.%A_%a.out

module load hmmer/3

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

N=$SLURM_ARRAY_TASK_ID
PEPDIR=pep
MARKERS=HMM/$HMM/markers_3.hmmb
CUTOFF=1e-10
OUT=search/$HMM
LIST=list # this is the list file

# can pass which file to process on cmdline too, eg bash jobs/01_hmmsearch.sh 1
if [ ! $N ]; then
  N=$1
fi

if [ ! $N ]; then
 echo "need to have a job id"
 exit;
fi

# number of processors to use set by PBS - can change or set this to a variable
# in config perhaps ?
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=1
fi

mkdir -p $OUT
# translate the number to a line number in the list of proteomes file to determine what to run
G=`sed -n ${N}p $LIST`

# convention is they all end in .aa.fasta - change this if not or make a variable
NM=`basename $G .aa.fasta`
echo "g=$G"

if [[ ! -f "$OUT/$NM.domtbl" || $PEPDIR/$G -nt $OUT/$NM.domtbl ]]; then
 hmmsearch --cpu $CPU -E $CUTOFF --domtblout $OUT/$NM.domtbl $MARKERS $PEPDIR/$G >& $OUT/$NM.log
else
 echo "skipping $NM - has already run"
fi
