#!/bin/bash

# could all this be replaced with snakemake/makefile??

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=hmmalign
#SBATCH --time=2:00:00
#SBATCH --output=logs/hmmaln.%A_%a.out

ALNCHUNK=1
OUTPEPEXT=aa.fa
OUTCDSEXT=cds.fa
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

while getopts c:f:i:s:H: OPT; do
  case $OPT in
    f) FORCE=${OPTARG} ;;
    i) IN=${OPTARG} ;;
    s) ALNCHUNK=${OPTARG} ;;
    H) HMM=${OPTARG} ;;
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

MARKERS=()
if [[ ${SLURM_ARRAY_TASK_ID} ]]; then
  # start when we have chunk ranges to run
  START=$(python -c "print( ((${SLURM_ARRAY_TASK_ID} - 1) * $ALNCHUNK) + 1)")
  END=$(python -c "print( ${SLURM_ARRAY_TASK_ID} * $ALNCHUNK )")
  sed -n "${START},${END}p" "$ALNFILES" | while read NM
  do
    MARKERS+=( $(basename "$IN" ".$OUTPEPEXT" ) )
  done
else
  MARKERS=( $(basename "$IN" ".$OUTPEPEXT") )
fi

MARKER="$(basename "$IN" ".$OUTPEPEXT")"
echo "IN=$IN gene=$MARKER"



for MARKER in "${MARKERS[@]}"
do
    echo "gene=$MARKER"
    echo " make -f $PHYLING_DIR/util/makefiles/Makefile.hmmalign SEQOVERLAP=$SEQOVERLAP RESOVERLAP=$RESOVERLAP HMM=$HMM $DIR/$MARKER.aa.clipkit"
    make -f $PHYLING_DIR/util/makefiles/Makefile.hmmalign SEQOVERLAP=$SEQOVERLAP RESOVERLAP=$RESOVERLAP HMMFOLDER=${DBDIR} HMM=$HMM $DIR/$MARKER.aa.clipkit 
    if [ -f $DIR/$MARKER.$OUTCDSEXT ]; then
	    make -f $PHYLING_DIR/util/makefiles/Makefile.hmmalign SEQOVERLAP=$SEQOVERLAP RESOVERLAP=$RESOVERLAP HMMFOLDER=${DBDIR} HMM=$HMM $DIR/$MARKER.cds.clipkit
    fi
done
