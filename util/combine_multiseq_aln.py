#!/usr/bin/env python3
import warnings
import os.path, argparse, sys, inspect, re, random

wrap = 60

script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,script_path+"/../")

import lib.libPHYling as PHYling

# this assumes input is in fasta format for simplicity

parser = argparse.ArgumentParser(description=help,add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-d','--dir', required=True,
                    help="The input directory for the alignments")
parser.add_argument('--ext', default="aa.clipkit",
                    help="Extension of the alignments to process (post trimming usually)")
parser.add_argument('--expected',default="expected_prefixes.lst",
                    help='File which lists the names of input prefixes expected, the union of all prefixes')
parser.add_argument('-o','--outfile','--out',required=True,
                    help='Output file for the concatenated alignments')
parser.add_argument('-v','--debug',action='store_true',
                    help='Debug messages printed')
parser.add_argument('--include',
                    help='Include file to specify a subset of input alignment files')
parser.add_argument('--seed',type=int,
                    help='Seed for random num generator for rand subselect')
parser.add_argument('--rand',type=int,
                    help='Select a random subset (this many) alignments') 

parser.add_argument('--moltype',default="PROT",
                    help='Molecule type listing the partition file (DNA or PROT)') 


parser.add_argument('-p','--partitions',
                    help='partitions file to save this result into') 

args = parser.parse_args(sys.argv[1:])

expected = []

partition_file = ""
if args.partitions:
    partition_file = args.partitions
else:
    baseout = os.path.splitext(args.outfile)
    partition_file = baseout[0] + ".partitions.txt"
    
with open(args.expected,"r") as fh:
    for line in fh:
        expected.append( re.sub(r'^>','',line.strip("\n")))

allowed = {}

if args.include:
    with open(args.include,"r") as fh:
        for line in fh:
            g = line.strip().split()
            allowed[g[0]] = 1

matrix = {}
partitions = []
last = 1
files = []

for file in os.listdir(args.dir):
    if file.endswith("."+args.ext):
        #if file == outfile:
        #    continue
        filearray = os.path.splitext(file)
        stem = filearray[0]
        if args.include:
            warnings.warn("checking if",stem,"is in the allowed set")
            if stem in allowed:
                files.append([stem,os.path.join(args.dir,file)])
        else:
            # pushing all files in since we did not set a restricted set
            # this is default behavior
            files.append([stem,os.path.join(args.dir,file)])

if args.rand:
    randfiles = random.shuffle(files)
    newfiles = []
    for n in range(args.rand):
        newfiles.append(randfiles[n])
    files = newfiles # does this do a proper copy?

# process every file, keep track of start and end for the
# partitions and make a concatenated alignment

for fname in files:
    stem     = fname[0]
    fullfile = fname[1]
    fstem    = re.sub(r'\.msa','',stem)
    seqs     = PHYling.read_fasta(fullfile)

    alnlen = len(seqs[0][1]) # first sequence, the alignseq is index [1]
    now = last + alnlen - 1
    
    partitions.append("%s, %s = %d-%d" % (args.moltype, fstem, last,now))
    last = now + 1

    seen = {}

    for seq in seqs:
        idstr = seq[0]
        seqstr = seq[1]
	
        if len(idstr) == 0:
            warnings.warn("no idstr for %s %s %s" % (idstr,fname, seqstr))
            exit()

        if len(seqstr) == 0:
            warnings.warn("no length for seq",idstr)
            exit()

        m = re.search(r'([^\|]+)\|',idstr)
        if m:
            idstr = m.group(1)

        # replace . to - in the alignment in case those are there
        seqstr = re.sub(r'\.','-',seqstr).upper()
        if idstr in matrix:
            matrix[idstr] += seqstr
        else:
            matrix[idstr] = seqstr

        seen[idstr] = 1

    for exp in expected:
        if exp not in seen:
            gapstr = '-' * alnlen
            if exp in matrix:
                matrix[exp] += gapstr
            else:
                matrix[exp] = gapstr

with open(args.outfile,"w") as fh:
    for seqid in matrix:
        data = matrix[seqid]

        fh.write(">%s\n"%(seqid))
        for i in range(0, len(data), wrap):
            fh.write(data[i:i + wrap] + "\n")

with open(partition_file,"w") as ofh:
    for p in partitions:
        ofh.write(p+"\n")
