#!/usr/bin/env bash

PEPDIR=pep
CDSDIR=cds
PEPEXT=aa.fasta
CDSEXT=cds.fasta
LOG_FOLDER=logs
if [ ! -f config.txt ]; then
 echo "Need a config.txt file"
 exit
fi

source ./config.txt
# make sure log folder is created
mkdir -p $LOG_FOLDER

if [ ! -d HMM/${HMM}  ]; then 
 echo "need a HMM defined in config.txt and this needs to be a folder in HMM/"
 echo "usually you checkout / download https://github.com/1KFG/Phylogenomics_HMMs/releases/latest and then "
 echo "symlink/move Phylogenomics_HMMs/HMM ./HMM"
 exit
fi

echo "make expected_prefixes.lst"
head -q -n1 $PEPDIR/*.$PEPEXT | awk -F\| '{print $1}' | awk '{print $1}' > expected_prefixes.lst

if [ ! -f prefix.tab ]; then
    echo "making prefix.tab"
    for file in $PEPDIR/*.${PEPEXT}
    do
	name=$(basename $file .${PEPEXT} | perl -p -e 's/\.\w+\.v\d+//;')
	pref=$(head -n1 $file | perl -p -e 's/^>([^\|\s]+).+/$1/')
	echo -e "$pref\t$name"
    done > prefix.tab
else
    echo "prefix.tab already exists, not updating"
fi
echo "make pepfile.lst for hmmsearch runs"
for f in $PEPDIR/*.$PEPEXT
do
    basename $f
done > pepfile.lst
