#!/usr/bin/bash -l
module load hmmer/3

for HMM in $(ls Kickxello_8/HMM3/*.hmm)
do
    name=$(basename $HMM .hmm)
    echo -n -e "$name\t"
    esl-seqstat -a pep/$name.fasta.pep | grep '^=' | awk '{print $3}' | ./make_lengths_from_seqs.pl
done > Kickxello_8/lengths_cutoff


for HMM in $(ls Kickxello_8/HMM3/*.hmm)
do
    name=$(basename $HMM .hmm)
    echo -n -e "$name\t"
    hmmsearch --domtbl $name.tbl $HMM $name.fasta.pep >& $name.hmmsearch
    grep -v '^#' $name.tbl | awk '{print $8}' | ./make_scores_cutoff.pl --adjust 0.80
done > Kickxello_8/scores_cutoff



