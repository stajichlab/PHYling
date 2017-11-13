#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=makeUnaln
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=makeunaln.%A.out

# DO NOT RUN WITH ARRAYJOBS

if [ $MODULESHOME ]; then
    module load cdbfasta
    module load hmmer
fi

PEPEXT=aa.fasta
CDSEXT=cds.fasta
ALN=aln
HMMSEARCH_OUTDIR=search

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

mkdir -p $ALN/$HMM

# rewrite with python?

# looks for .best files
# consider how this determines if it needs to run
# will not overwrite for now - how to force
perl $SCRIPTDIR/../util/make_unaln.pl -d $HMMSEARCH_OUTDIR/$HMM \
     -db pep -o $ALN/$HMM -ext $PEPEXT

if [ -d cds ]; then
    perl $SCRIPTDIR/../util/make_unaln.pl -d $HMMSEARCH_OUTDIR/$HMM \
	 -db cds -o $ALN/$HMM -ext $CDSEXT
fi
