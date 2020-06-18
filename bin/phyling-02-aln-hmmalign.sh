#!/bin/bash

# could all this be replaced with snakemake/makefile??

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=hmmalign
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmaln.%A_%a.out

OUTPEPEXT=aa.fa
RESOVERLAP=0.70
SEQOVERLAP=80
TRIMALSCHEME=-automated1
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Need a \"$CONFIG\" file"
    exit
fi

if [[ "$MODULESHOME" ]]; then
    module load hmmer/3
    module load trimal
    module load clipkit
    # module load BMGE # is this too slow?
fi

if [[ -z "$HMM_FOLDER" ]]; then
    echo "Need HMM_FOLDER set"
    exit
fi

if [[ ! -d "$HMM_FOLDER" ]]; then
    echo "HMM_FOLDER \"$HMM_FOLDER\" does not exist"
    exit
fi

ALN_OUTDIR="aln"

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

DIR="$ALN_OUTDIR/$HMM"
DBDIR="$HMM_FOLDER/$HMM/HMM3"

if [[ ${SLURM_ARRAY_TASK_ID} ]]; then
    IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ALNFILES")
fi

MARKER="$(basename "$IN" ".$OUTPEPEXT")"
echo "IN=$IN gene=$MARKER"

make -f $PHYLING_DIR/util/makefiles/Makefile.hmmalign SEQOVERLAP=$SEQOVERLAP RESOVERLAP=$RESOVERLAP $DIR/$MARKER.aa.clipkit

exit

#OUTFILE="$DIR/$MARKER.aa.msa"
#INFILE="$IN"

#if [[ "$FORCE" == "1" || ! -f "$OUTFILE" || "$IN" -nt "$OUTFILE" ]]; then
#    hmmalign --trim --amino -o "$OUTFILE" "$DBDIR/$MARKER.hmm" "$INFILE"
#fi

#INFILE="$OUTFILE" # last OUTFILE is new INFILE
#OUTFILE="$DIR/$MARKER.aa.clnaln"
#AA_ALN_NOTRIM="$OUTFILE"

#if [[ $FORCE == "1" || ! -f "$OUTFILE" || "$IN" -nt "$OUTFILE" ]]; then
#    TMP=$(mktemp)
#    esl-reformat --replace=x:- --gapsym=- -o "$TMP" afa "$INFILE"
#    perl -p -e 'if (! /^>/) { s/[ZBzbXx\*]/-/g }' "$TMP" > "$OUTFILE"
#    rm "$TMP"
#fi

#INFILE="$OUTFILE" # last OUTFILE is new INFILE
#OUTFILE="$DIR/$MARKER.aa.filter"

#if [[ $FORCE == "1" || ! -f $OUTFILE || $IN -nt $OUTFILE  ]]; then
#    trimal -resoverlap "$RESOVERLAP" \
#        -seqoverlap "$SEQOVERLAP" \
#        -in "$INFILE" \
#        -out "$OUTFILE"
#fi

#INFILE="$OUTFILE" # last OUTFILE is new INFILE
#OUTFILE="$DIR/$MARKER.aa.trim"

#if [[ $FORCE == "1" || ! -f "$OUTFILE" || "$IN" -nt "$OUTFILE" ]]; then
#    trimal "$TRIMALSCHEME" -fasta -in "$INFILE" -out "$OUTFILE"
#fi

#CDSFASTA="$DIR/$MARKER.cds.fa"
#CDSALN="$DIR/$MARKER.cdsaln"
#CDSTRIM="$DIR/$MARKER.cdsaln.trim"
#if [[ -f "$CDSFASTA" ]]; then
#    echo "processing CDS $CDSFASTA"
#    if [[ ! -f "$CDSALN" ]]; then
#        "$PHYLING_DIR/util/bp_mrtrans.pl" -if fasta -of fasta \
#                                -i "$AA_ALN_NOTRIM" \
#                                -s "$CDSFASTA" \
#                                -o "$CDSALN"
#        if [[ $BMGE ]]; then
#            java -jar "$BMGE" -t "CODON" -i "$CDSALN" -of "$CDSTRIM"
	#        else
#       rsync -a "$CDSALN" "$CDSTRIM"
#        fi
#    fi
#fi
