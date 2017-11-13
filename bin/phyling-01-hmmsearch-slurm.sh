#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=hmmsearch
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmsearch.%A.out

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

if [ ! $HMM ]; then
 echo "need to a config file to set the HMM folder name"
fi

N=$SLURM_ARRAY_TASK_ID
MARKERS=HMM/$HMM/markers_3.hmmb
OUT=$HMMSEARCH_OUT/$HMM

# can pass which file to process on cmdline too, eg bash jobs/01_hmmsearch.sh 1
if [ ! $N ]; then
  N=$1
fi

#if [ ! $N ]; then
# echo "need to have a job id"
# exit;
#fi

if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

mkdir -p $OUT
# translate the number to a line number in the list of proteomes file to determine what to run
#G=`sed -n ${N}p $LIST`
echo "$IN"
# convention is they all end in .aa.fasta - change this if not or make a variable
NM=`basename $IN.aa.fasta`
echo "g=$IN"

if [[ ! -f "$OUT/$NM.domtbl" || $PEPDIR/$IN -nt $OUT/$NM.domtbl ]]; then
 hmmsearch --cpu $CPU -E $HMMSEARCH_CUTOFF --domtblout $OUT/$NM.domtbl $MARKERS $PEPDIR/$IN >& $OUT/$NM.log
else
 echo "skipping $NM - has already run"
fi
