# PHYling tool
The unified PHYling pipeline for phylogenomic data collection from annotated genomes.

This is latest iteration of tool for using phylogenetically conserved markers to pull out informative 
gene or protein info from genomic and transcriptomic datasets in order to construct gene trees and species phylogenies.

The aligned markers can be extracted from protein sequences for phylogenetic analyses and also projected into coding sequence alignments for codon-based analyses for better resolution of recently diverged species. 

The assumptions in this approach are that the markers are generally single copy in genomes and taking best hit is sufficient first approximation for identifying orthologs. A separate file is parsed and file best_multihits which lists all the hits above the cutoff threshold for a given marker which can be used to assess duplication or attempt to incorporate paralogs into the analysis down the road. This approach can also be used for a quick-and-dirty quality check of your data by running a `wc -l search/*/*.best` and seeing a simple report of number of recovered markers.

Pay attention to how files are named and formatted in the pep folder. Currently expects a prefix embedded in the sequence ID e.g.
```
>Scer|YAL001C
```
If this prefix is missing the system cannot determine which taxon a sequence comes from when pulling data out of all the results in one go and to link these together for the concatenated result file. Check the `prefix.tab` file after running `PHYling init` to make sure it looks sensible with the prefix (not a whole gene name) assigned to each of the datafiles.

For example `prefix.tab` for the files in the [test directory](https://github.com/stajichlab/PHYling_unified/tree/main/test/pep) looks like (col1 is the prefix col2 is the species name gleaned from the filenanes, you can edit this file with full species names (includeing spaces if you want) to be used later in the phylogeny mapping from short prefixes back to full species names in the tree file:
```
Afum	Afum
B0A48	Rant
Scer	Scer
Spom	Spombe
Zymps1	Zymps1
```

The marker sets developed for this approach in fungi are available as part of the [1KFG Phylogenomics_HMMs](https://github.com/1KFG/Phylogenomics_HMMs) project resource and preferred use of the [BUSCO marker sets](https://busco-data.ezlab.org/v4/data/lineages/).

Usage
=====

You need in your path:
- `hmmsearch`, `hmmbuild` from [hmmer](http://hmmer.org/) version 3.x and the easel-tools part of the package (`esl-sfetch`)
- [clipkit](https://jlsteenwyk.com/ClipKIT/)
- python3
- to run multithreaded then install the [parallel](https://www.gnu.org/software/parallel/) unix tool `-q slurm` or if you are on a slurm supporting HPC can use the `-q slurm`
```
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
cd ..
curl -L -O ${BASEURL}/test/config.txt
# edit config.txt to indicate the prefix for this project and HMM marker set folder (eg fungi_odb10)
# set the PHYLING_DIR folder for now until we get this installable as a package
# set the number of processors for jobs and any job queue requirements for slurm in config.txt

# use the path to the installed PHYling or run this out of the checked-out folder - the folder where PHYling is installed
# needs to be specify as relative or absolute path in the config.txt
PHYling init
PHYling search -q parallel [or -q slurm]
PHYling aln -q parallel [or -q slurm]

```

Requirements
============
1. [HMMer](http://hmmer.org/) v3 - hmmsearch, hmmalign, esl-sfetch 
1. python v3
1. Perl (v5 with [BioPerl](http://bioperl.org) for tree rewriting if desired)
1. [clipkit](https://jlsteenwyk.com/ClipKIT/) for trimming
1. slurm for job queuing based job submission; [parallel](https://www.gnu.org/software/parallel/) for parallel analysis. Lacking these will simply run analysis in serial fashion and not able to take advantages of multiple cores/CPUs
1. [FastTree](http://www.microbesonline.org/fasttree/)/[IQTREE](http://www.iqtree.org/)/[RAxML](https://sco.h-its.org/exelixis/software.html) for downstream phylogenetic analyses
1. [muscle](https://www.drive5.com/muscle/) (optional) for de-novo multiple alignment instead of HMM guided hmmalign

Notes
==== 
1. You need to format your protein files with a prefix separated by a '|' e.g. -> see the test downloads above for working files
```
>Scer|YAL001W
```
1. the gene tree + ASTRAL coalescent is not integrated into the system yet, we run this separately - see this [example](https://github.com/stajichlab/Fusarium_Phylogenomics/blob/main/pipeline/06_pep_gene_trees.sh) which uses the Makefile approach to generate trees and this approach to [run ASTRAL](https://github.com/stajichlab/Fusarium_Phylogenomics/blob/main/pipeline/07_ASTRAL.sh)
2. the current strategy uses makefiles to solve dependency around pipeline steps
3. still need to add options to support additional parameters to the aligners and allow subset-filtering of markers for final concatenated file
3. The partitions file likely needs small changes for IQ-TREE or RAxML so it may not work directly depending on protein vs DNA trees
3. If you have coding sequences in the `cds` folder these will automatically be used to generate a back translated CDS alignment. The IDs in the protein file and CDS files have to be identical
4. training your own marker set is also possible but most busco sets are good starting place
