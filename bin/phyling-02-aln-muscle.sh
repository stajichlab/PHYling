#!/bin/bash

# could all this be replaced with snakemake/makefile??

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=muscle
#SBATCH --time=12:00:00
#SBATCH --output=logs/muscle.%A_%a.out

LOG_FOLDER=logs
QUEUEING=parallel
JOBCPU=1
TOTALCPU=4
OUTPEPEXT=aa.fa
OUTCDSEXT=cds.fa
#SCRIPTDIR=$(which PHYling | dirname)

if [ $MODULESHOME ]; then
    module load muscle
    module load trimal
#    module load java
#    module load BMGE # is this too slow?
fi

ALN_OUTDIR=aln
HMM_FOLDER=HMM

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM variable"
 exit
fi

#echo "$@"
while getopts c:f:i: option
do
 case "${option}"
 in
  c) CLEAN=${OPTARG};;
  f) FORCE=${OPTARG};;
  i) IN=${OPTARG};;
 esac
done

#echo "FORCE=$FORCE IN=$IN CLEAN=$CLEAN"

if [ ! $ALNFILES ]; then
    ALNFILES=alnlist.$HMM
fi

if [ ! -f $ALNFILES ]; then
    echo "expected an ALNFILES: $ALNFILES to exist"
    exit
fi

if [ ! $HMM ]; then
 echo "need to a config file to set the HMM folder name"
 exit
fi

DIR=${ALN_OUTDIR}/$HMM
DBDIR=${HMM_FOLDER}/${HMM}/HMM3

if [ ${SLURM_ARRAY_TASK_ID} ]; then
    IN=$(sed -n ${SLURM_ARRAY_TASK_ID}p $ALNFILES)
fi

if [ ${SLURM_CPUS_ON_NODE} ]; then
 CPU=${SLURM_CPUS_ON_NODE}
fi

marker=$(basename $IN .$OUTPEPEXT)
echo "IN=$IN gene=$marker"

OUTFILE=$DIR/$marker.aa.denovo.msa
INFILE=$IN
if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    muscle -in $IN -out $OUTFILE
fi

INFILE=$OUTFILE # last OUTFILE is new INFILE
OUTFILE=$DIR/$marker.aa.denovo.clean

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    esl-reformat --replace=\*:- --gapsym=- fasta $INFILE | esl-reformat --replace=x:- fasta > $OUTFILE
fi

INFILE=$OUTFILE # last OUTFILE is new INFILE
OUTFILE=$DIR/$marker.aa.denovo.trim

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    trimal -resoverlap 0.50 -seqoverlap 60 -in $INFILE -out $OUTFILE
    trimal -automated1 -fasta -in $INFILE -out $OUTFILE
fi

#if [ -f $DIR/$marker.cds.fasta ]; then
#    if [ ! -f $DIR/$marker.cdsaln.trim ]; then
#	$SCRIPTDIR/util/bp_mrtrans.pl -if clustalw -of fasta \
#				      -i $DIR/$marker.aln \
#				      -s $DIR/$marker.cds.fasta \
#				      -o $DIR/$marker.cdsaln
#	java -jar $BMGE -t CODON -i $DIR/$marker.cdsaln -of $DIR/$marker.cdsaln.trim
#    fi
#fi
