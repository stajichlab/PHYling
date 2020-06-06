#!/usr/bin/env python
import argparse, sys, csv
parser = argparse.ArgumentParser(description="Get top hit from HMM search",
                                 add_help=True)
parser.add_argument('-s','--scores',required=True,
                    help="HMM scores file cutoff (BUSCO scores_cutoff)")
parser.add_argument('-l','--lengths',required=True,
                    help="HMM lengths file cutoff (BUSCO lengths_cutoff)")
parser.add_argument('-i','--input',required=True,help="Input domtbl file")
parser.add_argument('-o','--output',required=False,help="Output multihit best report file")
args = parser.parse_args(sys.argv[1:])

#print(args)
#print(args.cutoff)
#print(args.input)
cutoffs = {}
with open(args.lengths,"r") as fh:
    for line in fh:
        (id,n,sd,mean) = line.strip().split("\t")
        if id in cutoffs:
            cutoffs[id]["length"] = [int(n),float(sd),float(mean)]
        else:
            cutoffs[id] = {"length": [int(n),float(sd),float(mean)]}

with open(args.scores,"r") as fh:
    for line in fh:
        (id,score) = line.strip().split("\t")
        #print(id,score)
        if id in cutoffs:
            cutoffs[id]["score"] = float(score)
        else:
            cutoffs[id] = {"score": score}


with open(args.input,"r") as fh:
    rowwriter = csv.writer(sys.stdout,delimiter="\t")
    seenbest = {}
    for line in fh:
        if line[0] != "#":
            line = line.strip("\n")
            row = line.split()
            t = row[0]
            tlen = row[2]
            q = row[3]
            evalue = row[6]
            score  = row[7]
            # we should add a component to check length
            if ( float(score) >= cutoffs[q]["score"] ):
                if q not in seenbest:
                    seenbest[q] = {t: ",".join([t,tlen, score, evalue])}
                elif t not in seenbest[q]:
                    seenbest[q][t] = ",".join([t,tlen, score, evalue])

    for s in seenbest:
        row = [s]
        for hit in seenbest[s]:
            row.append(seenbest[s][hit])
        rowwriter.writerow(row)
