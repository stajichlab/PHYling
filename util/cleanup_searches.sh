#!/bin/bash

source config.txt

for file in search/$HMM/*.log; do 
 b=$(basename $file .log); 
 if [ ! -f $PEPDIR/${b}.${PEPEXT} ]; then
  echo "rm search/$HMM/$b.*"
 fi
done
