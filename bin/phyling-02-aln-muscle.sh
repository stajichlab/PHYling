#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=muscle
#SBATCH --time=12:00:00
#SBATCH --output=logs/muscle.%A_%a.out

OUTPEPEXT=aa.fa
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Need a \"$CONFIG\" file"
    exit
fi

if [[ -z "$MODULESHOME" ]]; then
    module load muscle
    module load trimal
    #    module load java
    #    module load BMGE # is this too slow?
fi

while getopts c:f:i: OPT; do
    case $OPT in
          f) FORCE=${OPTARG};;
          i) IN=${OPTARG};;
    esac
done

[[ -z "$ALNFILES" ]] && ALNFILES="alnlist.$HMM"

if [[ ! -f "$ALNFILES" ]]; then
    echo "ALNFILES \"$ALNFILES\" does not exist"
    exit 1
fi

if [[ -z "$HMM" ]]; then
    echo "Need config file to set the HMM folder name"
    exit 1
fi

DIR="aln/$HMM"

if [[ ${SLURM_ARRAY_TASK_ID} ]]; then
    IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ALNFILES")
fi

MARKER=$(basename "$IN" ".$OUTPEPEXT")
echo "IN=$IN gene=$MARKER"

OUTFILE="$DIR/$MARKER.aa.denovo.msa"
INFILE="$IN"
if [[ "$FORCE" == "1" || ! -f "$OUTFILE" || "$IN" -nt "$OUTFILE" ]]; then
    muscle -in "$IN" -out "$OUTFILE"
fi

INFILE="$OUTFILE" # last OUTFILE is new INFILE
OUTFILE="$DIR/$MARKER.aa.denovo.clean"

if [[ $FORCE == "1" || ! -f "$OUTFILE" || "$IN" -nt "$OUTFILE" ]]; then
    esl-reformat --replace=\*:- --gapsym=- fasta "$INFILE" | \
        esl-reformat --replace=x:- fasta > "$OUTFILE"
fi

INFILE="$OUTFILE" # last OUTFILE is new INFILE
OUTFILE="$DIR/$MARKER.aa.denovo.trim"

if [[ "$FORCE" == "1" || ! -f "$OUTFILE" || "$IN" -nt "$OUTFILE" ]]; then
    trimal -resoverlap 0.50 -seqoverlap 60 -in "$INFILE" -out "$OUTFILE"
    trimal -automated1 -fasta -in "$INFILE" -out "$OUTFILE"
fi
