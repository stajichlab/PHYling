#!/usr/bin/env bash

#defaults which are overridden config.txt if set
QUEUEING=parallel
JOBCPU=1
ALN_OUTDIR=aln
TRIMALNEXT=trim
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Need a \"$CONFIG\" file"
    exit
fi
if [[ ! -d "$HMM_FOLDER" ]]; then
    echo "HMM_FOLDER \"$HMM_FOLDER\" does not exist"
    exit
fi

while getopts c:f:i: OPT; do
    case $OPT in
      m) METHOD=${OPTARG};;
      o) OUTDIR=${OPTARG};;
      q) QUEUEING=${OPTARG};;
    esac
done
