#!/usr/bin/env bash

#defaults which are overridden config.txt if set
QUEUEING=parallel
JOBCPU=1
TOTALCPU=4
OUTPEPEXT=aa.fa
ALN_OUTDIR=aln
ALNTOOL=hmmalign
EXPECTED=expected_prefixes.lst
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Need a \"$CONFIG\" file"
    exit
fi

JOBCPU=1 # force this to 1 since all of these steps are single threaded

if [[ -z "$HMM" ]]; then
    echo "Need config file to set the HMM folder name"
    exit 1
fi

[[ -z "$ALNFILES" ]] && ALNFILES="alnlist.$HMM"

while getopts t:c:f:q: OPT; do
    case $OPT in
         t) ALNTOOL=${OPTARG};;
         c) CLEAN=${OPTARG};;
         f) FORCE=${OPTARG};;
         q) QUEUEING=${OPTARG};;
    esac
done

if [[ ! -f "$ALNFILES" || "$CLEAN" == "1" || "$FORCE" == "1" ]]; then
    find "$ALN_OUTDIR/$HMM" -name "*.$OUTPEPEXT" > "$ALNFILES"
fi

SCRIPT_DIR=$(dirname "$0")
SUBJOB_SCRIPT="${SCRIPT_DIR}/phyling-02-aln-$ALNTOOL.sh"
COMBINE_SCRIPT="${SCRIPT_DIR}/phyling-03-aln-combine.sh"

if [[ $QUEUEING == "parallel" ]]; then
    JOBPARALLEL=$((TOTALCPU / JOBCPU))
    echo "Run parallel job $ALNTOOL"
    #echo "$JOBPARALLEL $SUBJOB_SCRIPT"
    parallel -j "$JOBPARALLEL" \
        "$SUBJOB_SCRIPT" \
        -f "$FORCE" \
        -c "$CLEAN" \
        -i {} < "$ALNFILES"

    if [[ $ALNTOOL == "muscle" ]]; then
	    "$COMBINE_SCRIPT" -x "$EXPECTED" --ext aa.denovo.trim
    else
	    "$COMBINE_SCRIPT" -x "$EXPECTED"
    fi

#    if [[ -f $CDSDIR/$ALLSEQNAME

elif [[ $QUEUEING == "slurm" ]]; then
    QUEUECMD=""
    [[ -n $QUEUE ]] && QUEUECMD="-p $QUEUE"

    ALNCT=$(wc -l "$ALNFILES" | awk '{print $1}')
    PHYLING_DIR=$(dirname "$0")

    echo "PHYLING_DIR is $PHYLING_DIR"
    ARRAY="1-$ALNCT"

    CMD="sbatch --ntasks $JOBCPU --nodes 1 $QUEUECMD --export=PHYLING_DIR=\"$PHYLING_DIR\" --export=FORCE=$FORCE --array=$ARRAY $SUBJOB_SCRIPT"

    SUBMIT_MSG=$(eval $CMD)
    SUBMIT_ID=$(echo "$SUBMIT_MSG" | awk '{print $4}')
    echo "SUBMIT_ID is $SUBMIT_ID; $SUBMIT_MSG"

    if [[ $ALNTOOL == "muscle" ]]; then
	echo "ready to run with $COMBINE_SCRIPT no extra ext"
        sbatch --depend=afterany:"$SUBMIT_ID" \
            $QUEUECMD \
            --export=EXT=aa.denovo.trim,EXPECTED=\"$EXPECTED\" \
            $COMBINE_SCRIPT
	
	if [ -f $CDSDIR/"cds_$ALLSEQNAME" ]; then
            sbatch --depend=afterany:"$SUBMIT_ID" \
		$QUEUECMD \
		--export=EXT=cds.denovo.trim,EXPECTED=\"$EXPECTED\" \
		"$COMBINE_SCRIPT"
	fi
    else
	echo "ready to run with $COMBINE_SCRIPT no extra ext"
        CMD="sbatch --depend=afterany:$SUBMIT_ID $QUEUECMD --export=EXPECTED=\"$EXPECTED\" $COMBINE_SCRIPT"
	eval $CMD

	if [ -f $CDSDIR/"cds_$ALLSEQNAME" ]; then
            CMD="sbatch --depend=afterany:$SUBMIT_ID $QUEUECMD --export=EXPECTED=$EXPECTED $COMBINE_SCRIPT"
	    eval $CMD
	fi
    fi
else
    echo "Run in serial"
    TMP=$(mktmp)
    find "$ALN_OUTDIR/$HMM" -name "*.$OUTPEPEXT" > "$TMP"
    while read -r FILE; do
        $SUBJOB_SCRIPT -f "$FORCE" -i "$FILE"
    done < "$TMP"
    rm "$TMP"
    # do combine here (how to avoid duplicate code with 
    # the parallel stuff above
fi
