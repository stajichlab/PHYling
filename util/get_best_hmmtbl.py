#!/usr/bin/env python
import argparse, sys

parser = argparse.ArgumentParser(description="Get top hit from HMM search",
                                 add_help=True)
parser.add_argument('-c','--cutoff',default='1e-20',type=float,
                    help="HMM parse cutoff")
parser.add_argument('-i','--input',help="Input domtbl file")
args = parser.parse_args(sys.argv[1:])

#print(args)
#print(args.cutoff)
#print(args.input)
with open(args.input,"r") as fh:
    seenbest = {}
    for line in fh:
        if line[0] != "#":
            line = line.strip("\n")
            row = line.split()
            t = row[0]
            q = row[3]
            evalue = row[6]
            if (float(evalue) <= args.cutoff and # better than cutoff
                (not q in seenbest)):    # and not already seen
                seenbest[q] = [q, t, evalue]
                
    for s in seenbest:
        print("\t".join(seenbest[s]))
