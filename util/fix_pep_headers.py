#!/usr/bin/env python3
"""
Author : kyclark
Date   : 2019-03-04
Purpose: Fix input peptide headers
"""

import argparse
import os
import sys
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """get command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Fix input peptide headers',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', metavar='str', help='Input file(s)', nargs='+')

    parser.add_argument(
        '-o',
        '--outdir',
        help='Output directory',
        metavar='str',
        type=str,
        required=True)

    return parser.parse_args()


# --------------------------------------------------
def warn(msg):
    """Print a message to STDERR"""
    print(msg, file=sys.stderr)


# --------------------------------------------------
def die(msg='Something bad happened'):
    """warn() and exit with error"""
    warn(msg)
    sys.exit(1)


# --------------------------------------------------
def main():
    """Make a jazz noise here"""
    args = get_args()
    out_dir = args.outdir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    for i, file in enumerate(args.file, start=1):
        basename = os.path.basename(file)
        print('{:3}: {}'.format(i, basename))

        out_file = os.path.join(out_dir, basename)
        out_fh = open(out_file, 'wt')
        prefix = basename.split('.')[0] # STR0027594.proteins.fa => STR0027594
        for rec in SeqIO.parse(file, 'fasta'):
            rec.id = '|'.join([prefix, rec.id])
            rec.description = ''
            SeqIO.write(rec, out_fh, 'fasta')

        out_fh.close()

    print('Done')


# --------------------------------------------------
if __name__ == '__main__':
    main()
