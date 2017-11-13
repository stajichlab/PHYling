#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=hmmalign
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --output=hmmalign.%A_%a.out
#SBATCH -p intel

if [ $MODULESHOME ]; then
    module load trimal
    module load hmmer/3
    module load java
#    module load BMGE # is this too slow?
fi

ALN_OUTDIR=aln
HMM_FOLDER=HMM
LIST=alnlist.$MARKER # this is the list file
SCRIPTDIR=$(dirname $0)

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set HMM variable"
 exit
fi

if [ ! $HMM ]; then
 echo "need to a config file to set the HMM folder name"
fi
DIR=${ALN_OUTDIR}/$HMM
DBDIR=${HMM_FOLDER}/${HMM}/HMM3

if [ ! $IN ]; then
  IN=$1
fi

marker=`basename $IN .aa.fasta`
echo "marker is $marker for gene $IN"

if [ ! -f $DIR/$marker.msa ]; then
    hmmalign --trim --amino $DBDIR/$marker.hmm $IN > $DIR/$marker.msa
    # hmmalign --trim --amino $DBDIR/$marker.hmm $DIR/$marker.aa.fasta | perl -p -e 's/^>(\d+)\|/>/' > $DIR/$marker.msa
fi

if [ ! -f $DIR/$marker.aln ]; then
    esl-reformat --replace=\*:-  --gapsym=- clustal $DIR/$marker.msa > $DIR/$marker.1.aln
    esl-reformat --replace=x:- clustal $DIR/$marker.1.aln > $DIR/$marker.aln
fi

if [ ! -f $DIR/$marker.msa.trim ]; then
    trimal -resoverlap 0.50 -seqoverlap 60 -in $DIR/$marker.aln -out $DIR/$marker.msa.filter
    trimal -automated1 -fasta -in $DIR/$marker.msa.filter -out $DIR/$marker.msa.trim 
fi

if [ -f $DIR/$marker.cds.fasta ]; then
    if [ ! -f $DIR/$marker.cdsaln.trim ]; then
	$SCRIPTDIR/util/bp_mrtrans.pl -if clustalw -of fasta \
				      -i $DIR/$marker.aln \
				      -s $DIR/$marker.cds.fasta \
				      -o $DIR/$marker.cdsaln
#	java -jar $BMGE -t CODON -i $DIR/$marker.cdsaln -of $DIR/$marker.cdsaln.trim
    fi
fi
