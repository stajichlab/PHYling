#!/bin/bash

# could all this be replaced with snakemake/makefile??

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=hmmalign
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmaln.%A_%a.out

LOG_FOLDER=logs
QUEUEING=parallel
JOBCPU=1
TOTALCPU=4
OUTPEPEXT=aa.fa
OUTCDSEXT=cds.fa
RESOVERLAP=0.50
SEQOVERLAP=60
TRIMALSCHEME=-automated1

#SCRIPTDIR=$(which PHYling | dirname)

if [ $MODULESHOME ]; then
    module load hmmer/3
    module load trimal
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

OUTFILE=$DIR/$marker.aa.msa
INFILE=$IN

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    hmmalign --trim --amino -o $OUTFILE $DBDIR/$marker.hmm $INFILE 
fi

INFILE=$OUTFILE # last OUTFILE is new INFILE
OUTFILE=$DIR/$marker.aa.clnaln

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    esl-reformat --replace=x:- --gapsym=- -o $OUTFILE.tmp afa $INFILE
    perl -p -e 'if (! /^>/) { s/[ZBzbXx\*]/-/g }' $OUTFILE.tmp > $OUTFILE
    rm $OUTFILE.tmp
fi

INFILE=$OUTFILE # last OUTFILE is new INFILE
OUTFILE=$DIR/$marker.aa.filter

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    trimal -resoverlap 0.50 -seqoverlap 60 -in $INFILE -out $OUTFILE
fi

INFILE=$OUTFILE # last OUTFILE is new INFILE
OUTFILE=$DIR/$marker.aa.filter

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    trimal -resoverlap $RESOVERLAP -seqoverlap $SEQOVERLAP -in $INFILE -out $OUTFILE
fi

INFILE=$OUTFILE # last OUTFILE is new INFILE
OUTFILE=$DIR/$marker.aa.trim

if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
    trimal $TRIMALSCHEME -fasta -in $INFILE -out $OUTFILE
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
