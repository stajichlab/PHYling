# PHYling_unified
Unified PHYling pipeline for phylogenomic data collection from annotated genomes.

This is latest iteration of tool for using phylogenetically conserved markers to pull out informative 
gene or protein info from genomic and transcriptomic datasets in order to construct gene trees and species phylogenies.

The aligned markers can be extracted from protein sequences for phylogenetic analyses and also projected into coding sequence alignments for codon-based analyses for better resolution of recently diverged species. 

The assumptions in this approach are that the markers are generally single copy in genomes and taking best hit is sufficient first approximation for identifying orthologs. 

The marker sets developed for this approach in fungi are available as part of the [1KFG Phylogenomics_HMMs](https://github.com/1KFG/Phylogenomics_HMMs) project resource.

Usage
=====


Requirements
============
1. [HMMer](http://hmmer.org/) v3 - hmmsearch, hmmalign, esl-sfetch 
1. python v3
1. Perl (v5 with [BioPerl](http://bioperl.org) for tree rewriting if desired)
1. slurm for job queuing based job submission; [parallel](https://www.gnu.org/software/parallel/) for parallel analysis. Lacking these will simply run analysis in serial fashion and not able to take advantages of multiple cores/CPUs
1. [FastTree](http://www.microbesonline.org/fasttree/)/[IQTREE](http://www.iqtree.org/)/[RAxML](https://sco.h-its.org/exelixis/software.html) for downstream phylogenetic analyses
1. [muscle](https://www.drive5.com/muscle/) (optional) for de-novo multiple alignment instead of HMM guided hmmalign
