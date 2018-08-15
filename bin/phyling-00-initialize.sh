#!/usr/bin/env bash
#SBATCH --time 2:00:00 --out=init.log -J init
QUERYDBS=pepfile.lst
PEPDIR=pep
CDSDIR=cds
INPEPEXT=aa.fasta
LOG_FOLDER=logs
EXPECTED=expected_prefixes.lst
if [ ! -f config.txt ]; then
 echo "Need a config.txt file"
 exit
fi

mkdir -p $LOG_FOLDER

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
head -q -n1 $PEPDIR/*.$INPEPEXT | awk -F\| '{print $1}' | awk '{print $1}' > $EXPECTED

if [ ! -f prefix.tab ]; then
    echo "making prefix.tab"
    for file in $PEPDIR/*.${INPEPEXT}
    do
	name=$(basename $file .${INPEPEXT} | perl -p -e 's/\.\w+\.v\d+//;')
	pref=$(head -n1 $file | perl -p -e 's/^>([^\|\s]+).+/$1/')
	echo -e "$pref\t$name"
    done > prefix.tab
else
    echo "prefix.tab already exists, not updating"
fi
echo "make $QUERYDBS for hmmsearch runs"
for f in $PEPDIR/*.$INPEPEXT
do
    basename $f
done > $QUERYDBS
