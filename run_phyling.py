#!/usr/bin/env python3
"""Run PHYling"""

import argparse
import configparser
import inspect
import json
import lib.libPHYling as PHYling
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
import zipfile
from urllib.request import urlopen


# --------------------------------------------------
def get_args():
    """get args"""
    parser = argparse.ArgumentParser(
        description='Run PHYling',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'command',
        metavar='COMMAND',
        help='Command',
        default='')

    parser.add_argument(
        '-c',
        '--config',
        help='Config file name',
        metavar='str',
        type=str,
        default='config.txt')

    sub_parsers = parser.add_subparsers(help='Command help')
    sub_parsers.add_parser('init', help='Setup data directory, needs HMM dir')

    dl_parser = sub_parsers.add_parser(
        'download', help='Download HMM for pre-defined phylogenomic markers')
    dl_parser.add_argument(
        '-t', help='HMM name', metavar='str', type=str, default='fungi')

    sub_parsers.add_parser(
        'search', help='Search HMM set to the genomes in pep')

    sub_parsers.add_parser(
        'aln', help='Construct unaligned FASTA of protein/CDS, align')

    sub_parsers.add_parser(
        'superaln', help='Concatenate gene alignments into a superalignment')

    phylo_parser = sub_parsers.add_parser(
        'phylo', help='Run Phylogenetic reconstruction on superalignment')

    phylo_parser.add_argument(
        '-t',
        help='raxml, iqtree, fasttree',
        metavar='str',
        type=str,
        default='raxml')

    genetree_parser = sub_parsers.add_parser(
        'genetrees',
        help='Run phylogenetic reconstruction on gene from marker')

    genetree_parser.add_argument(
        '-t', help='raxml, fasttree', metavar='str', type=str, default='raxml')

    genetree_parser.add_argument(
        '-m', help='aa, cds', metavar='str', type=str, default='aa')

    coalesce_parser = sub_parsers.add_parser(
        'coalesce', help='Run ASTRAL gene tree reconciliation software')

    coalesce_parser.add_argument(
        '-t',
        help='raxml, iqtree, fasttree',
        metavar='str',
        type=str,
        default='raxml')

    coalesce_parser.add_argument(
        '-m', help='aa, cds', metavar='str', type=str, default='aa')

    coalesce_parser.add_argument(
        '-b',
        help='Sample bootstrap files if generated (RAXML)',
        action='store_true')

    sub_parsers.add_parser('citation', help='Print citation')

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
    """main"""
    version = '0.2'  #phyling release version
    logging.basicConfig()

    args = get_args()
    subprog = args.command
    print('>>>>>>> {}'.format(subprog))

    script_path = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe())))

    # this code loads up some config files which list where the HMM libraries
    # can be downloaded from
    URL_file = os.path.join(script_path, "lib", "urls.json")
    Messages_file = os.path.join(
        script_path, "lib",
        "messages.json")  # abstract this later based on different languages
    HMMs_URL = {}
    Messages = {}
    with open(URL_file, "r") as jsonfile:
        jsonToPython = json.loads(jsonfile.read())
        HMMs_URL = jsonToPython['HMMs']

    with open(Messages_file, "r") as jsonfile:
        Messages = json.loads(jsonfile.read())['en']

    default_help = Messages['commands']['default'] % version

    if subprog == "version":
        print("PHYling v%s" % version)
        sys.exit(0)
    elif subprog == "citation":
        print(Messages["citation"])
        sys.exit(0)

    config = get_config(args.config)

    if subprog == 'init' or subprog == 'initialize' or subprog == 'setup':
        parser = args.init_parser
        cmd = os.path.join(script_path, 'bin', 'phyling-00-initialize.sh')
        subprocess.call(cmd)
    elif subprog == 'download':
        help = Messages['commands']['download'] % (sys.argv[1], version)
        parser = argparse.ArgumentParser(
            description=help,
            add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)

        parser.add_argument('-t', '--type', default='fungi')
        args = parser.parse_args(arguments)
        url = ""
        outfile = args.type + "_HMMs.zip"

        #HEREHERE
        if args.type in HMMs_URL:
            url = HMMs_URL[args.type]

            if not os.path.exists(outfile):
                with urlopen(url) as response, open(outfile, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
                    unzippath = args.type + "_HMMs"

                if not os.path.exists(unzippath):
                    zip_ref = zipfile.ZipFile(outfile, 'r')
                    zip_ref.extractall(unzippath)

                if not os.path.exists(config["HMM_FOLDER"]):
                    os.makedirs(config["HMM_FOLDER"])

                for subdir in os.listdir(unzippath):
                    dlong = os.path.join(unzippath, subdir)
                    sublist = os.listdir(dlong)
                    if 'HMM' in sublist:
                        for hmmfolder in os.listdir(
                                os.path.join(dlong, 'HMM')):
                            fromf = os.path.join(dlong, 'HMM', hmmfolder)
                            print(fromf)
                            shutil.move(fromf, config["HMM_FOLDER"])
                            print(
                                "%s HMMs installed. URL=%s Outfile=%s Dest=%s"
                                % (args.type, url, outfile,
                                   config["HMM_FOLDER"]))

    elif re.match("(hmm|search)", subprog):
        help = Messages['commands']['search'] % ('search', version)
        parser = argparse.ArgumentParser(
            description=help,
            add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument(
            '-f',
            '--force',
            action='store_true',
            help="force overwriting of files when running")
        parser.add_argument(
            '-q', '--queueing', help="Queueing parallel, serial, or slurm")
        args = parser.parse_args(arguments)
        if args.queueing is not None:
            config["QUEUEING"] = args.queueing

        if args.force:
            print(args.force)
            searchfolder = os.path.join(config["HMMSEARCH_OUTDIR"],
                                        config["HMM"])
            for file in os.listdir(searchfolder):
                if (file.endswith(".domtbl") or file.endswith(".log")
                        or file.endswith("." + config["BESTHITEXT"])):
                    os.unlink(os.path.join(searchfolder, file))

        cmd = os.path.join(script_path, 'bin', 'phyling-01-hmmsearch.sh')
        subprocess.call([
            cmd,
            "-f",
            "%d" % int(args.force),
            "-q",
            "%s" % (config["QUEUEING"]),
        ])

    elif re.match("aln", subprog):
        help = Messages['commands']['aln'] % ('aln', version)
        parser = argparse.ArgumentParser(
            description=help,
            add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        # hmmer or muscle for multiple alignment?
        parser.add_argument(
            '-t', '--type', action='store_true', default="hmmalign")

        # force clean DB and align
        parser.add_argument('-f', '--force', action='store_true')

        # clean align folder (remake starting files)
        parser.add_argument('-c', '--cleanaln', action='store_true')
        # override queue option from config file
        parser.add_argument(
            '-q', '--queueing', help="Queueing parallel, serial, or slurm")

        args = parser.parse_args(arguments)

        if args.queueing is not None:
            config["QUEUEING"] = args.queueing

        outdir = config["PEPDIR"]
        if "TEMP" in config:
            outdir = os.path.join(config["TEMP"], config["PREFIX"])
            if not os.path.isdir(outdir):
                os.makedirs(outdir)

        pep_db = os.path.join(outdir, config["ALLSEQNAME"])

        PHYling.init_allseqdb(config["PEPDIR"], pep_db, config["INPEPEXT"],
                              config["INDEXING"], args.force)

        searchdir = os.path.join(config["HMMSEARCH_OUTDIR"], config["HMM"])
        alndir = os.path.join(config["ALN_OUTDIR"], config["HMM"])

        # parse the best hit files, make ortholog table and write out genes
        # to a single file per ortholog  - first do the proteins
        #print("make unaln called with",searchdir,pep_db)
        PHYling.make_unaln_files(searchdir, config["BESTHITEXT"],
                                 config['HMMSEARCH_CUTOFF'], pep_db, alndir,
                                 config["OUTPEPEXT"], args.cleanaln,
                                 config["INDEXING"], int(config["TOTALCPU"]))
        # now re-parse the best hit files, make ortholog table and write out genes
        # to a single file per ortholog for the coding sequence files
        # assuming there is a CDS folder (skip if not)

        if os.path.exists(config["CDSDIR"]):
            outdir = config["CDSDIR"]
            if "TEMP" in config:
                outdir = os.path.join(config["TEMP"], config["PREFIX"])
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)

            cds_db = os.path.join(outdir, "cds_" + config["ALLSEQNAME"])
            PHYling.init_allseqdb(config["CDSDIR"], cds_db, config["INCDSEXT"],
                                  config["INDEXING"], args.force)

            PHYling.make_unaln_files(searchdir, config["BESTHITEXT"],
                                     config['HMMSEARCH_CUTOFF'], cds_db,
                                     alndir, config["OUTCDSEXT"],
                                     args.cleanaln, config["INDEXING"],
                                     int(config["TOTALCPU"]))

        # either force or cleanaln flag sufficient to regenerate alignment files
        # do we sub-divide by type (muscle,hmmalign here)
        cmd = os.path.join(script_path, 'bin', 'phyling-02-aln.sh')
        #print(cmd)
        subprocess.call([
            cmd,
            "-t",
            args.type,
            "-c",
            "%d" % (int(args.cleanaln)),
            "-f",
            "%d" % (int(args.force)),
            "-q",
            "%s" % (config["QUEUEING"]),
        ])

    elif subprog == "phylo":
        help = Messages['commands']['phylo'] % (sys.argv[1], version)
        parser = argparse.ArgumentParser(
            description=hqelp,
            add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        # hmmer or muscle for multiple alignment?
        args = parser.parse_args(arguments)

    elif re.match(r"genetree", subprog):
        help = Messages['commands']['genetrees'] % (sys.argv[1], version)
        parser = argparse.ArgumentParser(
            description=help,
            add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        # hmmer or muscle for multiple alignment?
        args = parser.parse_args(arguments)
    elif re.match(r"coal", subprog) or subprog == "astral":
        help = Messages['commands']['coalesce'] % (sys.argv[1], version)
        parser = argparse.ArgumentParser(
            description=help,
            add_help=True,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        args = parser.parse_args(arguments)

    else:
        print("%s option not recognized" % sys.argv[1])
        print(default_help)


# --------------------------------------------------
def get_config(config_file):
    config = {
        'HMM_FOLDER': 'HMM',
        'HMM': 'AFTOL70',
        'PEPDIR': 'pep',
        'INPEPEXT': 'aa.fasta',  # input data
        'OUTPEPEXT': 'aa.fa',  # output alignment files
        'CDSDIR': 'cds',
        'INCDSEXT': 'cds.fasta',  # input data
        'OUTCDSEXT': 'cds.fa',  # output alignment files
        'BESTHITEXT': 'best',
        'HMMSEARCH_CUTOFF': 1e-30,
        'HMMSEARCH_OUTDIR': 'search',
        'ALN_OUTDIR': 'aln',
        'ALLSEQNAME': 'allseq',
        'LANGUAGE': 'en',
        'PREFIX': 'PHY',
        'QUEUEING': 'parallel',
        'INDEXING': 'sfetch',
    }

    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                line = re.sub(r"\s*#.+$", "", line)
                keyval = line.split("=")
                config[keyval[0]] = keyval[1]
    else:
        die('Bad --config "{}" file'.format(config_file))

    return config


# --------------------------------------------------
if __name__ == '__main__':
    main()
