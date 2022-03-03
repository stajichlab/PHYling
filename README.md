# PHYling tool
The unified PHYling pipeline for phylogenomic data collection from annotated genomes.

This is latest iteration of tool for using phylogenetically conserved markers to pull out informative 
gene or protein info from genomic and transcriptomic datasets in order to construct gene trees and species phylogenies.

The aligned markers can be extracted from protein sequences for phylogenetic analyses and also projected into coding sequence alignments for codon-based analyses for better resolution of recently diverged species. 

The assumptions in this approach are that the markers are generally single copy in genomes and taking best hit is sufficient first approximation for identifying orthologs. 

The marker sets developed for this approach in fungi are available as part of the [1KFG Phylogenomics_HMMs](https://github.com/1KFG/Phylogenomics_HMMs) project resource and preferred use of the [BUSCO marker sets](https://busco-data.ezlab.org/v4/data/lineages/).

Usage
=====

```
module load hmmer/3
module unload miniconda
module load miniconda3

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
BASEURL=https://raw.githubusercontent.com/stajichlab/PHYling_unified/main/
mkdir -p pep
cd pep
for SP in Afum Rant Zymps1 Scer Spom
do
  curl -L -O ${BASEURL}/test/pep/${SP}.aa.fasta
done

# edit config.txt to indicate the prefix for this project and HMM marker set folder (eg fungi_odb10)
# set the PHYLING_DIR folder for now until we get this installable as a package
curl -L -O ${BASEURL}/test/config.txt

PHYling init
PHYling search
PHYling aln

```

Requirements
============
1. [HMMer](http://hmmer.org/) v3 - hmmsearch, hmmalign, esl-sfetch 
1. python v3
1. Perl (v5 with [BioPerl](http://bioperl.org) for tree rewriting if desired)
1. slurm for job queuing based job submission; [parallel](https://www.gnu.org/software/parallel/) for parallel analysis. Lacking these will simply run analysis in serial fashion and not able to take advantages of multiple cores/CPUs
1. [FastTree](http://www.microbesonline.org/fasttree/)/[IQTREE](http://www.iqtree.org/)/[RAxML](https://sco.h-its.org/exelixis/software.html) for downstream phylogenetic analyses
1. [muscle](https://www.drive5.com/muscle/) (optional) for de-novo multiple alignment instead of HMM guided hmmalign
