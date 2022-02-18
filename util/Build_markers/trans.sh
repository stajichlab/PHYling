#!?usr/bin/bash -l

module unload miniconda2
module load miniconda3
module load muscle
module load hmmer/3

perl trans.pl
mkdir -p Kickxello_8/HMM3
for file in $(ls pep/*.pep)
do
	base=$(basename $file .fasta.pep)
	muscle < $file > $file.aln
	hmmbuild -n $base Kickxello_8/HMM3/$base.hmm $file.aln
done
