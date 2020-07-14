#!/usr/bin/env python3

import argparse, sys, csv, os

parser = argparse.ArgumentParser(description="Integrate prefix.tab with additional names (replacing those prefix.tab)",
                                 add_help=True)
parser.add_argument('-i','--intable',required=True,type=argparse.FileType('r'),
                    help="tab delimited 2 cols: Prefix<tab>LongName")
parser.add_argument('--dir',required=False,default="search",help="Top folder searches are stored in")
parser.add_argument('--hmm',required=False,default="fungi_odb10",help="Folder HMM Name to search for input")
args = parser.parse_args(sys.argv[1:])

starting = csv.reader(args.intable,delimiter="\t")
prefixes = {}

for row in starting:
    pref = row[0]
    name = row[1]
    if pref in prefixes:
        print("found %s in row %s, but already stored for %s, prefix input file is not unique?"%(pref,name,prefixes[pref][0]))

    prefixes[pref] = [name,0]

indir=os.path.join(args.dir,args.hmm)
for fname in os.listdir(indir):
    if fname.endswith(".best"):
        with open(os.path.join(indir,fname),"r") as fh:
            besthits = csv.reader(fh,delimiter="\t")
            count = 0
            prefix = ""
            for row in besthits:
                if count == 0:
                    name = row[1]
                    prefix = name.split('|')[0]
                count += 1
            if not prefix in prefixes:
                print("could not find prefix %s in parsed table of names (%s)"%(prefix,fname))
                prefixes[prefix] = [prefix,0]
            prefixes[prefix][1] = count

for prefix in sorted(prefixes.keys(), key=lambda x: prefixes[x][0]):
    if prefixes[prefix][1] > 0:
        print("%-20s\t%s\t%d"%(prefix,prefixes[prefix][0],prefixes[prefix][1]))
