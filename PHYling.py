#!/usr/bin/env python3

import sys, os, subprocess, inspect, tarfile, shutil, argparse
import urllib.request, configparser, re, distutils.util
import json, logging
import zipfile

from urllib.request import urlopen

#from pathlib import Path
#phylingpath = os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))),"lib")
#print(phylingpath)

import lib.libPHYling as PHYling

logging.basicConfig()

# TODO: general variables and paths need to be set and initialized
# here at the beginning

# user specific config file - the nane of this should be overrideable
Config_file = 'config.txt'

version = '0.2.6'  #phyling release version

# sys.argv parsing should occur to get/set an alternative config file name

script_path = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, script_path)

# some hard coded variables

# set some defaults which can be overridden by config.txt
# currently assume that input
config = {
    'HMM_FOLDER': 'HMM',
    'HMM': 'AFTOL70',
    'PEPDIR': 'pep',
    'INPEPEXT': 'aa.fasta',  # input data
    'OUTPEPEXT': 'aa.fa',  # output alignment files
    'CDSDIR': 'cds',
    'INCDSEXT': 'cds.fasta',  # input data
    'OUTCDSEXT': 'cds.fa',  # output alignment files
    'BESTHITEXT': 'best_multi',
    'INCLUDE_MULTI': 'False',
    'HMMSEARCH_CUTOFF': 1e-20,
    'HMMSEARCH_OUTDIR': 'search',
    'ALN_OUTDIR': 'aln',
    'ALLSEQNAME': 'allseq',
    'LANGUAGE': 'en',
    'PREFIX': 'PHY',
    'QUEUEING': 'parallel',
    'INDEXING': 'sfetch',
}

if os.path.exists(Config_file):
    with open(Config_file, "r") as configparse:
        for line in configparse:
            line = line.strip()
            if line.startswith("#"):
                continue
            line = re.sub(r"\s*#.+$", "", line)
            keyval = line.split("=")
            config[keyval[0]] = keyval[1]
else:
    print(
        "Expecting a config.txt file to define the input folders, outgroup, and HMM marker set"
    )

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
    Messages = json.loads(jsonfile.read())[config["LANGUAGE"]]

default_help = Messages['commands']['default'] % version

# if no arguments, print out help message
if len(sys.argv) <= 1:
    print(default_help)
    exit(1)

# first argument in the sysargs is going to be the subprogram name
# e.g. init, download, ...
subprog = sys.argv[1].lower()
arguments = sys.argv[2:]

if subprog == 'init' or subprog == 'initialize' or subprog == 'setup':
    help = Messages['commands']['initialize'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(
        description=help,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', description="Extra debugging")
    # parse the rest of the argument
    args = parser.parse_args(arguments)
    cmd = os.path.join(script_path, 'bin', 'phyling-00-initialize.sh')
    subprocess.call(cmd)
elif subprog == 'download':
    help = Messages['commands']['download'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(
        description=help,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', description="Extra debugging")
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
                    for hmmfolder in os.listdir(os.path.join(dlong, 'HMM')):
                        fromf = os.path.join(dlong, 'HMM', hmmfolder)
                        print(fromf)
                        shutil.move(fromf, config["HMM_FOLDER"])
                        print("%s HMMs installed. URL=%s Outfile=%s Dest=%s" %
                              (args.type, url, outfile, config["HMM_FOLDER"]))

elif re.match("(hmm|search)", subprog):
    help = Messages['commands']['search'] % ('search', version)
    parser = argparse.ArgumentParser(
        description=help,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help="Extra debugging")
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
        searchfolder = os.path.join(config["HMMSEARCH_OUTDIR"], config["HMM"])
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
    parser.add_argument('-v', '--verbose', action='store_true',help="Extra debugging")
    parser.add_argument(
        '-t', '--type', action='store_true', default="hmmalign",help="What aligner - hmmsearch or muscle")

    # force clean DB and align
    parser.add_argument('-f', '--force', action='store_true', help="Force re-running of steps even if output file exists")
    parser.add_argument('-m', '--multi', action='store_true',
                        default=distutils.util.strtobool(config["INCLUDE_MULTI"]),
                        help="Keep multi-hits?")

    # clean align folder (remake starting files)
    parser.add_argument('-c', '--cleanaln', action='store_true',help="Clean alignments before running")
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
    #print("make unaln called with",searchdir,pep_db,args.multi)

    PHYling.make_unaln_files(
        searchdir, config["BESTHITEXT"], config['HMMSEARCH_CUTOFF'], pep_db,
        alndir, config["OUTPEPEXT"], args.cleanaln, config["INDEXING"],
        int(config["TOTALCPU"]),args.multi)
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
                                 config['HMMSEARCH_CUTOFF'], cds_db, alndir,
                                 config["OUTCDSEXT"], args.cleanaln,
                                 config["INDEXING"], int(config["TOTALCPU"],
                                 args.multi))

    # either force or cleanaln flag sufficient to regenerate alignment files
    # do we sub-divide by type (muscle,hmmalign here)

    cmd = os.path.join(script_path, 'bin', 'phyling-02-aln.sh')
    print("calling aln step: {} -t {} -c {} -q {}".format(cmd,args.type,int(args.cleanaln),config["QUEUEING"]))

    subprocess.call([
        cmd,
        "-t",
        args.type,
        "-c",
        "%d" % (int(args.cleanaln)),
        "-q",
        "%s" % (config["QUEUEING"]),
    ])

elif subprog == "phylo":
    help = Messages['commands']['phylo'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(
        description=help,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # hmmer or muscle for multiple alignment?
    args = parser.parse_args(arguments)

elif re.match(r"genetree", subprog):
    help = Messages['commands']['genetree'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(
        description=help,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # iqtree or raxml
    parser.add_argument('-q', '--queueing', help="Queueing parallel, serial, or slurm")

    parser.add_argument('-m', '--method', default="fasttree",
                        help="Tree Building method")

    parser.add_argument('-o', '--outputdir', default="gene_trees",
                        help="Output Gene Tree Folder")

    args = parser.parse_args(arguments)
    cmd = os.path.join(script_path, 'bin', 'phyling-05-genetree.sh')
    print(cmd)
    subprocess.call([
        cmd,
        "-m",
        args.method,
        "-o",
        args.outputdir,
        "-q",
        "%s" % (config["QUEUEING"]),
    ])

elif re.match(r"coal", subprog) or subprog == "astral":
    help = Messages['commands']['coalesce'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(
        description=help,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args(arguments)

elif subprog == "version":
    print("PHYling v%s" % version)
elif subprog == "citation":
    print(Messages["citation"])
else:
    print("%s option not recognized" % sys.argv[1])
    print(default_help)
