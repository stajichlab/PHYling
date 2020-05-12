#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=hmmsearch
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmsearch.%A_%a.out

[[ $MODULESHOME ]] && module load hmmer/3

QUERYDBS=pepfile.lst
PEPDIR=pep
PEPEXT=aa.fasta
HMM_FOLDER=HMM
HMMSEARCH_CUTOFF=1e-10
HMMSEARCH_OUT=search
CPU=2
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Expected config file \"$CONFIG\" to set HMM variable"
    exit 1
fi

if [[ ! $PHYLING_DIR ]]; then
    echo "Need a PHYLING_DIR defined in \"$CONFIG\$ or env variable"
    exit 1
fi

if [[ ! $HMM ]]; then
    echo "Need config file or env to set the HMM folder name"
    exit 1
fi

MARKERS="$HMM_FOLDER/$HMM/markers_3.hmmb"
OUT="$HMMSEARCH_OUT/$HMM"
[[ ! -d "$OUT" ]] && mkdir -p "$OUT"

if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$QUERYDBS")
elif [[ ! $IN ]]; then
    # can pass which file to process on cmdline too,
    # eg bash jobs/01_hmmsearch.sh 1
    IN=$1
fi

[[ ${SLURM_CPUS_ON_NODE} ]] && CPU=${SLURM_CPUS_ON_NODE}

IN=$(basename "$IN")
NM=$(basename "$IN" .$PEPEXT)
echo "g=$IN NM=$NM"

INFILE="$PEPDIR/$IN"
OUTFILE1="$OUT/$NM.domtbl"
OUTFILE2="$OUT/$NM.log"

if [[ $FORCE == "1" || ! -f $OUTFILE1 || $INFILE -nt $OUTFILE1  ]]; then
    hmmsearch --cpu "$CPU" \
        -E "$HMMSEARCH_CUTOFF" \
        --domtblout "$OUTFILE1" \
        "$MARKERS" \
        "$INFILE" >& "$OUTFILE2"

else
    echo "Skipping $NM - has already run"
fi

OUTFILEBEST="$OUT/$NM.best"
if [[ "$FORCE" == "1" || ! -f "$OUTFILEBEST" || "$OUTFILE1" -nt "$OUTFILEBEST" ]]; then
    "${PHYLING_DIR}/util/get_best_hmmtbl.py" \
        -c $HMMSEARCH_CUTOFF \
        --input "$OUTFILE1" > "$OUTFILEBEST"
fi

OUTFILEBESTMULTI="$OUT/$NM.best_multi"
if [[ "$FORCE" == "1" || ! -f "$OUTFILEBESTMULTI" || "$OUTFILE1" -nt "$OUTFILEBESTMULTI" ]]; then
    "${PHYLING_DIR}/util/get_best_by_cutoffs_hmmtbl.py" \
        -s $HMM_FOLDER/$HMM/scores_cutoff -l $HMM_FOLDER/$HMM/lengths_cutoff \
        --input "$OUTFILE1" > "$OUTFILEBESTMULTI"
fi
