#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 16gb 
module load hmmer/3
module unload miniconda3
module load clipkit

DIR=HMM
MARKER=fungi_odb10
URL=https://busco-data.ezlab.org/v4/data/lineages/fungi_odb10.2020-09-10.tar.gz

if [ ! -d $DIR/fungi_odb10 ]; then
    mkdir -p $DIR/$MARKER
    curl -C - -O $URL
    FILE=$(basename $URL)
    tar zxf $FILE
    
    mv $MARKER/scores_cutoff $MARKER/lengths_cutoff $DIR/$MARKER
    mv $MARKER/hmms $DIR/$MARKER/HMM3
    cat $DIR/$MARKER/HMM3/*.hmm > $DIR/$MARKER/markers_3.hmm
    hmmconvert -b $DIR/$MARKER/markers_3.hmm > $DIR/$MARKER/markers_3.hmmb
    
    rm -rf fungi_odb10.2020-09-10.tar.gz $MARKER
fi

../PHYling init
../PHYling search -q parallel
../PHYling aln -q parallel
