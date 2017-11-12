#!/usr/bin/bash
perl -i -p -e 's/>jgi\|([^\|]+)\|(\d+)\|(\S+)/>$1|$1_$2 $3/' $1
