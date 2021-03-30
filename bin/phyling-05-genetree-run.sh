#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=2G
#SBATCH --job-name=genetree
#SBATCH --time=12:00:00
#SBATCH --output=logs/gene_tree.%A_%a.out

#defaults which are overridden config.txt if set
ALN_OUTDIR=aln
TRIMALNEXT=trim
CONFIG="config.txt"

if [[ -f "$CONFIG" ]]; then
    source "$CONFIG"
else
    echo "Need a \"$CONFIG\" file"
    exit
fi

if [[ -z "$MODULESHOME" ]]; then
    module load RAxML/8.2.12
    module load fasttree/2.1.11
    module load IQ-TREE/2.1.1
fi

if [[ ! -d "$HMM_FOLDER" ]]; then
    echo "HMM_FOLDER \"$HMM_FOLDER\" does not exist"
    exit
fi

while getopts c:f:i: OPT; do
    case $OPT in
      i) INPUT=${OPTARG};;
      o) OUTDIR=${OPTARG};;
      f) FORCE=${OPTARG};;
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

if [[ ${SLURM_ARRAY_TASK_ID} ]]; then
    IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ALNFILES")
else
  echo "IN not set becuse 'SLURM_ARRAY_TASK_ID' not set"
  exit
fi

MARKER="$(basename "$IN" ".$OUTPEPEXT.$TRIMALNEXT")"
echo "IN=$IN gene=$MARKER"

INFILE="$DIR/$MARKER.aa.$TRIMALNEXT"
OUTFILE="$MARKER"

make -f $PHYLING_DIR/util/makefiles/Makefile.trees $DIR/$MARKER.aa.clipkit.FT.tre
make -f $PHYLING_DIR/util/makefiles/Makefile.trees $DIR/$MARKER.cds.clipkit.FT.tre
