#!/usr/bin/bash -l

#defaults which are overridden config.txt if set
QUEUEING=parallel
JOBCPU=1
TOTALCPU=4
ALNCHUNK=1
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

while getopts t:c:s:q: OPT; do
  case $OPT in
    t) ALNTOOL=${OPTARG} ;;
    c) CLEAN=${OPTARG} ;;
    s) ALNCHUNK=${OPTARG} ;;
    q) QUEUEING=${OPTARG} ;;
  esac
done

if [[ ! -f "$ALNFILES" || "$CLEAN" == "1" ]]; then
  find "$ALN_OUTDIR/$HMM" -name "*.$OUTPEPEXT" > "$ALNFILES"
fi

SCRIPT_DIR=$(dirname "$0")
SUBJOB_SCRIPT="${SCRIPT_DIR}/phyling-02-aln-$ALNTOOL.sh"
COMBINE_SCRIPT="${SCRIPT_DIR}/phyling-03-aln-combine.sh"

if [[ $QUEUEING == "parallel" ]]; then
  JOBPARALLEL=$((TOTALCPU / JOBCPU))
  cat $ALNFILES | parallel -j "$JOBPARALLEL" \
    "$SUBJOB_SCRIPT" \
    -c "$CLEAN" \
    -i {} 

  if [[ $ALNTOOL == "muscle" ]]; then
	  echo "This may not be well supported - most testing has focused on HMMALIGN"
    "$COMBINE_SCRIPT" -x "$EXPECTED" -e aa.denovo.trim -t aa
    if [ -d $CDSDIR ]; then
      "$COMBINE_SCRIPT" -x "$EXPECTED" -e cds.denovo.trim -t cds
    fi
  else
    "$COMBINE_SCRIPT" -x "$EXPECTED" -e "aa.clipkit" -t aa
    if [ -d $CDSDIR ]; then
      "$COMBINE_SCRIPT" -x "$EXPECTED" -e "cds.clipkit" -m "DNA" -t cds
    fi
  fi


elif [[ $QUEUEING == "slurm" ]]; then
  QUEUECMD=""
  [[ -n $QUEUE ]] && QUEUECMD="-p $QUEUE"

  ALNCT=$(wc -l "$ALNFILES" | awk '{print $1}')
  ALNCTDIV=$(python -c "print(int($ALNCT / $ALNCHUNK) + 1)")
  PHYLING_DIR=$(dirname "$0")
  #echo "DEBUG: ALNCT is $ALNCT ALNCTDIV is $ALNCTDIV and ALNFILES was $ALNFILES"
  echo "DEBUG: PHYLING_DIR is $PHYLING_DIR"
  ARRAY="1-$ALNCTDIV"
  CMD="sbatch --ntasks $JOBCPU --nodes 1 $QUEUECMD --export=PHYLING_DIR=\"$PHYLING_DIR\" --array=$ARRAY $SUBJOB_SCRIPT"
  echo "DEBUG: $CMD"
  SUBMIT_MSG=$(eval $CMD)
  SUBMIT_ID=$(echo "$SUBMIT_MSG" | awk '{print $4}')
  echo "SUBMIT_ID is $SUBMIT_ID; $SUBMIT_MSG"

  if [[ $ALNTOOL == "muscle" ]]; then
    echo "DEBUG: eady to run with $COMBINE_SCRIPT no extra ext"
    sbatch --depend=afterany:"$SUBMIT_ID" \
      $QUEUECMD \
      --export=EXT=aa.denovo.clipkit,TYPE=aa,EXPECTED=\"$EXPECTED\" \
      $COMBINE_SCRIPT

    if [ -d $CDSDIR ]; then
      sbatch --depend=afterany:"$SUBMIT_ID" \
        $QUEUECMD \
        --export=EXT=cds.denovo.clipkit,TYPE=CDS,EXPECTED=\"$EXPECTED\" \
        "$COMBINE_SCRIPT"
    fi
  else
    echo "ready to run with $COMBINE_SCRIPT no extra ext"
    CMD="sbatch --depend=afterany:$SUBMIT_ID $QUEUECMD --export=EXT=aa.clipkit,TYPE=aa,EXPECTED=\"$EXPECTED\" $COMBINE_SCRIPT"
    eval $CMD

    if [ -d $CDSDIR ]; then
      CMD="sbatch --depend=afterany:$SUBMIT_ID $QUEUECMD --export=EXT=cds.clipkit,MOLTYPE=DNA,EXPECTED=\"$EXPECTED\",TYPE=CDS $COMBINE_SCRIPT"
      eval $CMD
    fi
  fi
else
  echo "Run in serial"
  TMP=$(mktmp)
  find "$ALN_OUTDIR/$HMM" -name "*.$OUTPEPEXT" > "$TMP"
  while read -r FILE; do
    $SUBJOB_SCRIPT -i "$FILE" -c "$CLEAN"
  done < "$TMP"
  rm "$TMP"
  "$COMBINE_SCRIPT" -x "$EXPECTED" -e "aa.clipkit" -t aa
  "$COMBINE_SCRIPT" -x "$EXPECTED" -e "cds.clipkit" -m DNA -t cds
fi
