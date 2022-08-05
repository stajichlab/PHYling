#!/usr/bin/bash -l

#SBATCH --time 2:00:00 --out=init.log -J init

QUERYDBS=pepfile.lst
PEPDIR=pep
INPEPEXT=aa.fasta
LOG_FOLDER=logs
EXPECTED=expected_prefixes.lst
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Need a \"$CONFIG\" file"
    exit
fi

mkdir -p "$LOG_FOLDER"


# make sure log folder is created
[[ ! -d "$LOG_FOLDER" ]] && mkdir -p "$LOG_FOLDER"

HMM_DIR="$HMM_FOLDER/$HMM"

if [[ ! -d "$HMM_DIR" ]]; then 
    echo "Bad HMM_DIR \"$HMM_DIR\""
    echo "need a HMM defined in config.txt and this needs to be a folder in HMM/"
    echo "usually you checkout / download https://github.com/1KFG/Phylogenomics_HMMs/releases/latest and then "
    echo "symlink/move Phylogenomics_HMMs/HMM ./HMM"
    exit
fi

echo "make expected_prefixes.lst"
head -q -n1 $PEPDIR/*.$INPEPEXT | awk -F\| '{print $1}' | awk '{print $1}' > $EXPECTED

if [[ ! -f prefix.tab ]]; then
    echo "Making prefix.tab"
    for FILE in $PEPDIR/*.${INPEPEXT}; do
        NAME=$(basename "$FILE" ".$INPEPEXT" | perl -p -e 's/\.\w+\.v\d+//;')
        PREF=$(head -n1 "$FILE" | perl -p -e 's/^>([^\|\s]+).+/$1/')
        echo -e "$PREF\t$NAME"
    done > prefix.tab
else
    echo "prefix.tab already exists, not updating"
fi

echo "Make $QUERYDBS for hmmsearch runs"

for FILE in $PEPDIR/*.$INPEPEXT; do
    basename "$FILE"
done > "$QUERYDBS"
