#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=hmmsearch
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmsearch.%A_%a.out

if [ $MODULESHOME ]; then
    module load hmmer/3
fi
QUERYDBS=pepfile.lst
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
 exit
fi

MARKERS=${HMM_FOLDER}/$HMM/markers_3.hmmb
OUT=$HMMSEARCH_OUT/$HMM
mkdir -p $OUT

if [ ${SLURM_ARRAY_TASK_ID} ]; then
    IN=$(sed -n ${SLURM_ARRAY_TASK_ID}p $QUERYDBS)
elif [ ! $IN ]; then
# can pass which file to process on cmdline too, eg bash jobs/01_hmmsearch.sh 1
  IN=$1
fi

if [ ${SLURM_CPUS_ON_NODE} ]; then
 CPU=${SLURM_CPUS_ON_NODE}
fi

IN=$(basename $IN)
NM=$(basename $IN .$PEPEXT)
echo "g=$IN NM=$NM"

INFILE=$PEPDIR/$IN
OUTFILE1=$OUT/$NM.domtbl
OUTFILE2=$OUT/$NM.log

if [[ $FORCE == "1" || ! -f $OUTFILE1 || $INFILE -nt $OUTFILE1  ]]; then

 hmmsearch --cpu $CPU -E $HMMSEARCH_CUTOFF --domtblout $OUTFILE1 $MARKERS $INFILE >& $OUTFILE2

else
 echo "skipping $NM - has already run"
fi

OUTFILEBEST=$OUT/$NM.best
if [[ $FORCE == "1" || ! -f $OUTFILEBEST || $OUTFILE1 -nt $OUTFILEBEST ]]; then
    ${PHYLING_DIR}/util/get_best_hmmtbl.py -c $HMMSEARCH_CUTOFF --input $OUTFILE1 > $OUTFILEBEST
fi
