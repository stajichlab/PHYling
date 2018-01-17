#!/usr/bin/env python3

import os.path, argparse
import lib.libPHYling as PHYling

# this assumes input is in fasta format for simplicity

parser = argparse.ArgumentParser(description=help,add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-d','--dir', action='store_true',required=true,
                    help="The input directory for the alignments")
parser.add_argument('--ext', default="msa.trim",
                    help="Extension of the alignments to process (post trimming usually)")
parser.add_argument('--expected',default="expected_prefixes.lst",
                    help='File which lists the names of input prefixes expected, the union of all prefixes')
parser.add_argument('-o','--outfile','--out',required=true,
                    help='Output file for the concatenated alignments')
parser.add_argument('-v','--debug',action='store_true',
                    help='Debug messages printed')
parser.add_argument('--include',
                    help='Include file to specify a subset of input alignment files')
parser.add_argument('--seed',type=int,
                    help='Seed for random num generator for rand subselect')
parser.add_argument('--rand',type=int,
                    help='Select a random subset (this many) alignments') 

args = parser.parse_args(sys.argv)

for file in os.listdir(args.dir):
    if file.endswith("."+args.ext):
        print(file)
        
#PHYling.read_fasta


