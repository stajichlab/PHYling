#!/usr/bin/env python3

import os
import re
import subprocess
import logging

from multiprocessing.dummy import Pool as ThreadPool

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# todo
# - need to deal with how to save databases for different BUSCO lineages in same file


def get_CDS_with_exonerate(lineage, genomedb, proteinfile, tmpfolder="tmp", exonerateapp="exonerate", window=50):
    #  index the genome database in fasta format for quick lookup and retrieval

    genomename = os.path.basename(genomedb)
    species = re.sub(r'\.dna\.fasta', '', genomename)
    buscoId = re.sub(r'\.faa', '', os.path.basename(proteinfile))
    idxfile = os.path.join(tmpfolder, "{}.idx".format(genomename))
    idx = SeqIO.index_db(idxfile, genomedb, "fasta")
    pepseq = 0
    for seq in SeqIO.parse(proteinfile, "fasta"):
        pepseq = seq
        break
    if not pepseq:
        print("no protein in file {}".format(proteinfile))
        return
    chromname, location = pepseq.id.split(":")
    start, end = location.split("-")
    start = int(start)
    end = int(end)
    if chromname not in idx:
        print("no {} in genome db {}".format(chromname, genomedb))
        return
    # make sure start < end
    if start > end:
        start, end = end, start
    start -= window
    # make sure we are within bounds
    if start < 1:
        start = 1

    end += window

    if end > len(idx[chromname]):
        end = len(idx[chromname])
    chromloc = "{}_{}_{}".format(chromname, start, end)
    chromfile = os.path.join(tmpfolder, "{}.dna".format(chromloc))

    if not os.path.exists(chromfile):
        subseq = SeqRecord(idx[chromname][start:end].seq,
                           id=chromloc, description="")
        SeqIO.write(subseq, chromfile, "fasta")

    exonerate_command = "{} -m p2g --bestn 1 --ryo '>%ti %tcb_%tce\n%tcs' --showalignment false --showvulgar false --refine full -q {} -t {}".format(
        exonerateapp, proteinfile, chromfile)
    # print(exonerate_command)
    log_exonerate = ""
    with subprocess.Popen(exonerate_command, shell=True, stdout=subprocess.PIPE) as exonerateproc:
        inseq = 0
        CDSstr = ""
        chromloc = ""
        for linein in exonerateproc.communicate():
            if linein is None:
                break
            formattedline = linein.decode('utf-8')
            log_exonerate = formattedline
            for line in formattedline.splitlines():
                if line.startswith(">"):
                    inseq = 1
                    m = re.match(r'^>(\S+)\_(\d+)_(\d+)', line)
                    if m:
                        seqID     = m.group(1)
                        seq_start = m.group(2)
                        seq_end   = m.group(3)
                        if seq_end < seq_start:
                            seq_end,seq_start = seq_start,seq_end
                        chromloc = '{}_{}_{}'.format(seqID, seq_start, seq_end)
                    else:
                        print("cannot match seqid from {}".format(line))
                        break
                elif line.startswith("-- completed"):
                    inseq = 0
                elif inseq:
                    CDSstr += line


        if CDSstr:
            CDS = SeqRecord(Seq(CDSstr), id=buscoId, description=chromloc)
            pep = SeqRecord(CDS.translate().seq, id=buscoId, description=chromloc)

            for peprecord in SeqIO.parse(proteinfile, "fasta"):
                if abs(len(peprecord) - len(pep)) > 2:
                    print("Warning: Inferred translated CDS length ({}) for {} {} {} is not same length as input pep ({})".format(len(pep),
                            buscoId, species, chromloc, len(peprecord)))
                    print(">{}\n{}".format(buscoId,CDSstr))
                    print("\t{}".format(log_exonerate))
                    print("==\n")

            return {'lineage': lineage,
                    'busco':   buscoId,
                    'species': species,
                    'CDS':     CDS,
                    'pep':     pep}
        else:
            return {'lineage': lineage,
                    'busco':   buscoId,
                    'species': species}


# Main program here

FoldersFromBusco = ["single_copy_busco_sequences"]

parser = argparse.ArgumentParser(
    description='Reorganize BUSCO single-copy genes into MultiFasta for PHYling')

parser.add_argument('-i', '--BUSCO', default='BUSCO',
                    help='BUSCO folder')

parser.add_argument('-o', '--aln', default='aln',
                    help='alignment output folder')

parser.add_argument('--pepdir', default='pep',
                    help='protein directory')
parser.add_argument('--pepdbextension', default='aa.fasta',
                    help='pep db file extension')


parser.add_argument('--cdsdir', default='cds',
                    help='protein directory')
parser.add_argument('--cdsdbextension', default='cds.fasta',
                    help='CDS db file extension')

parser.add_argument('--genomedir', default='genomes',
                    help='genomes directory')

parser.add_argument('--genomeextension', default='dna.fasta',
                    help='genome file extension')

parser.add_argument('--pepextension', default='aa.fa',
                    help='alned file extension')

parser.add_argument('--cdsextension', default='cds.fa',
                    help='alned file extension')


parser.add_argument('--exonerate', default='exonerate',
                    help='exonerate tool')

parser.add_argument('--temp', default='tmp',
                    help='temp folder')

parser.add_argument('--sfetch', default='esl-sfetch',
                    help='esl-sfetch path')

parser.add_argument('-p', '--prefixfile', default='prefix.tab',
                    help='prefixes')

parser.add_argument('-t', '--threads', default='4', type=int,
                    help='Number of parallel threads to use')

parser.add_argument('-a', '--arrayjob', required=False, type=int,
                    help='If running array jobs, only run this specific Nth BUSCO folder')
parser.add_argument('--debug', const=True, default=False, nargs='?',
                    help='Debug steps')
parser.add_argument('--force', const=True, default=False, nargs='?',
                    help='Force writing output DB if already exists')

parser.add_argument('--fragmented', const=True, default=False, nargs='?',
                    help='Include fragemented sequences')
args = parser.parse_args()

if args.fragmented:
    FoldersFromBusco.appedn("fragmented_busco_sequences")

# make directories for storing alignments and protein files
for d in [args.pepdir, args.cdsdir, args.aln, args.temp]:
    if not os.path.exists(d):
        os.mkdir(d)

species_seqs = {}
orthologs = {}
sp_folders = []

pool = ThreadPool(args.threads)
BUSCOdir = sorted(os.listdir(args.BUSCO))

if args.arrayjob is not None and args.arrayjob > 0:
    dirsToRun = [BUSCOdir[args.arrayjob - 1]]
else:
    dirsToRun = BUSCOdir

for spdir in dirsToRun:
    spdirpath = os.path.join(args.BUSCO, spdir)
    pepdb = os.path.join(args.pepdir, "{}.{}".format(
        spdir, args.pepdbextension))
    cdsdb = os.path.join(args.cdsdir, "{}.{}".format(
        spdir, args.cdsdbextension))
    if ((os.path.exists(cdsdb) and os.path.getsize(cdsdb) > 0) or
            (os.path.exists(pepdb) and os.path.getsize(pepdb) > 0)):
        if args.force:
            print(
                "{} or {} already exist but overwriting because of --force".format(cdsdb, pepdb))
        else:
            print("skipping because {} or {} already exists and not forcing".format(
                cdsdb, pepdb))
            continue

    sp_folders.append(spdir)
    idxfile = os.path.join(args.temp, "{}.dna.fasta.idx".format(spdir))
    genome_file = os.path.join(
        args.genomedir, "{}.{}".format(spdir, args.genomeextension))
    idx = SeqIO.index_db(idxfile, genome_file, "fasta")

    if spdir not in species_seqs:
        species_seqs[spdir] = {'CDS': [],
                               'pep': []}
    run_CDS_pep_gather = []
    for buscorun in os.listdir(spdirpath):
        if buscorun.startswith("run_"):
            for folderName in FoldersFromBusco:
                seqFolder = os.path.join(
                    spdirpath, buscorun, "busco_sequences", folderName)
                BUSCO_lineage = buscorun.replace("run_", "")

                for orthseqfile in os.listdir(seqFolder):
                    if not orthseqfile.endswith(".faa"):
                        continue
                    orthologname = orthseqfile.replace(".faa", "")
                    if orthologname not in orthologs:
                        orthologs[orthologname] = {BUSCO_lineage: {}}
                    orthologfile = os.path.join(seqFolder, orthseqfile)
#                    print("{} -> {}".format(genome_file,orthologfile))
                    run_CDS_pep_gather.append(
                        [BUSCO_lineage, genome_file, orthologfile, args.temp, args.exonerate])

    print("There are {} processes to run".format(len(run_CDS_pep_gather)))
    results = None
    if args.debug:
        results = pool.starmap(get_CDS_with_exonerate, run_CDS_pep_gather[0:5])
    else:
        results = pool.starmap(get_CDS_with_exonerate, run_CDS_pep_gather)

    for res in results:  # list

        name = res["busco"]
        buscolineage = res["lineage"]
        species = res["species"]

        if 'CDS' not in res:
            print("No CDS returned for {} {}".format(name, species))
            continue

        species_seqs[species]["CDS"].append(
            SeqRecord(res["CDS"].seq, id=name, description=res["CDS"].description))
        species_seqs[species]["pep"].append(
            SeqRecord(res["pep"].seq, id=name, description=res["pep"].description))

        orthologs[name][buscolineage][species] = {"CDS": species_seqs[species]["CDS"][-1],
                                                  'pep': species_seqs[species]["pep"][-1]}
pool.close()
pool.join()

prefixes = {}
prefixes_rev = {}
if os.path.exists(args.prefixfile):
    with open(args.prefixfile, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("Prefix"):
                row = line.split("\t")
                prefixes[row[0]] = row[1]
                prefixes_rev[row[1]] = row[0]

for sp in sp_folders:
    #       print("spname is {}".format(sp))
    prefix = ""
    longname = sp
    longname = re.sub(r'_$', '', longname)
    spmod = re.sub(r'(cf|sp)\._', '', longname)
    name = spmod.split("_", 2)
    name = spmod.split("_", 2)
    if len(name) == 1:
        prefix = name[0]
    elif len(name) ==2 or len(name[2]) == 0:
        prefix = name[0][0] + name[0][1] + name[1][0:4]
    elif len(name) > 2:
        prefix =  name[0][0] + name[0][1] + name[1][0:4] + name[2]
    else:
        prefix = spmod
    prefix = re.sub(r'-','_',prefix)
    print(prefix,longname)
    if prefix not in prefixes:
        prefixes[prefix] = sp
        prefixes_rev[sp] = prefix

    with open(args.prefixfile, "wt") as ofh:
        for p in sorted(prefixes):
            ofh.write("\t".join([p, prefixes[p]]) + "\n")

for sp in species_seqs:
    seqsrename = []
    pepdb = os.path.join(args.pepdir, "{}.{}".format(sp, args.pepdbextension))
    cdsdb = os.path.join(args.cdsdir, "{}.{}".format(sp, args.cdsdbextension))
    for seqrec in species_seqs[sp]["pep"]:
        seqsrename.append(SeqRecord(seqrec.seq, id="{}|{}".format(
            prefixes_rev[sp], seqrec.id), description=seqrec.description))
    SeqIO.write(seqsrename, pepdb, "fasta")

    seqsrename = []
    for seqrec in species_seqs[sp]["CDS"]:
        seqsrename.append(SeqRecord(seqrec.seq, id="{}|{}".format(
            prefixes_rev[sp], seqrec.id), description=seqrec.description))

    SeqIO.write(seqsrename, cdsdb, "fasta")
