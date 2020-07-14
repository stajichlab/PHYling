#!/usr/bin/env bash

QUERYDBS=pepfile.lst
JOBCPU=2
TOTALCPU=4
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Expected config file \"$CONFIG\" to set HMM variable"
    exit 1
fi

SCRIPT_DIR=$(dirname "$0")
SUBJOB_SCRIPT="${SCRIPT_DIR}/phyling-01-hmmsearch-run.sh"
PHYLING_DIR=$(dirname $SCRIPT_DIR)

while getopts f:q: OPT; do
    case $OPT in
	f) FORCE=${OPTARG};;
        q) QUEUEING=${OPTARG};;
    esac
done

if [[ $QUEUEING == "parallel" ]]; then
    JOBPARALLEL=$((TOTALCPU / JOBCPU))
    echo "Run parallel job hmmsearch"
    parallel -j "$JOBPARALLEL" "$SUBJOB_SCRIPT" < "$QUERYDBS"
elif [[ $QUEUEING == "slurm" ]]; then
    QUEUECMD=""
    [[ $QUEUE ]] && QUEUECMD="-p $QUEUE"

    FILECOUNT=$(wc -l "$QUERYDBS" | awk '{print $1}')
    ARRAY="1-$FILECOUNT"
    echo "sbatch --ntasks=$JOBCPU --nodes=1 $QUEUECMD --export=PHYLING_DIR="$PHYLING_DIR" --array=$ARRAY $SUBJOB_SCRIPT"
    sbatch --ntasks=$JOBCPU \
        --nodes=1 \
        $QUEUECMD \
        --export=PHYLING_DIR="$PHYLING_DIR" \
	       --array=$ARRAY \
        "$SUBJOB_SCRIPT"
else
    echo "Run in serial"
    while read -r FILE; do
        $SUBJOB_SCRIPT "$FILE"
    done < "$QUERYDBS"
fi
