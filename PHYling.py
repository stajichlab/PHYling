#!/usr/bin/env python3

#wrapper for PHYling

import sys, os, subprocess, inspect, tarfile, shutil, argparse, urllib.request, configparser, re
from pathlib import Path

import logging

Config_file = 'config.txt'
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,script_path)

HMMs_URL = { 'fungi': 'https://github.com/1KFG/Phylogenomics_HMMs/releases/download/v1.0/fungal_HMMs-1.0.zip',
}

config = {}
if os.path.exists(Config_file):
    with open(Config_file,"r") as configparse:
        for line in configparse:
            line = line.strip()
            line = re.sub(r"\s*#.+$","",line)
            keyval = line.split("=")
            config[keyval[0]] = keyval[1]
    if not('HMM_FOLDER' in config):
        config["HMM_FOLDER"] = 'HMM'    
else:
    print("Expecting a config.txt file to define the input folders, outgroup, and HMM marker set")
    sys.exit(1)
    

version = '0.1'
default_help = """
Usage:       PHYling <command> <arguments>
version:     %s

Description: PHYling is a package of scripts to extract phylogenomic markers
             Dependencies:  HMMER3 MUSCLE RAxML

init:        Setup or recreate initialize files
download:    Download HMM markers

Written by Jason Stajich (2014-2017) jason.stajich[at]ucr.edu or jasonstajich.phd[AT]gmail.com
 Started in https://github.com/1KFG/Phylogenomics and https://github.com/stajichlab/phyling
        """ % version

if len(sys.argv) > 1:
    if sys.argv[1].lower() == 'init' or sys.argv[1].lower() == 'initialize' or sys.argv[1].lower() == 'setup':
        help = """
Usage:       PHYling %s <arguments>
version:     %s

Description: Script will setup initial data directory and files. Expects an HMM folder
             config.txt and 'pep' folder as specified in config.txt file. An optional 'cds' folder
""" % (sys.argv[1], version)
        arguments = sys.argv[2:]
        parser = argparse.ArgumentParser(description=help,add_help=True,formatter_class=argparse.RawDescriptionHelpFormatter)
        args = parser.parse_args(arguments)
        cmd = os.path.join(script_path, 'bin', 'phyling-00-initialize.sh')
        subprocess.call(cmd)

    if sys.argv[1].lower() == 'download':
        help = """
Usage:       PHYling %s <arguments>
version:     %s

Description: Download HMM files for pre-defined phylogenomic markers. 

Arguments:   -t fungi [default] 

""" % (sys.argv[1], version)
        arguments = sys.argv[2:]
        parser = argparse.ArgumentParser(description=help,add_help=True,formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-t','--type',default='fungi')
        args = parser.parse_args(arguments)
        url = ""
        outfile = args.type + "_HMMs.zip"
        if args.type in HMMs_URL:
            url = HMMs_URL[args.type]

            if not os.path.exists(outfile):
                with urllib.request.urlopen(url) as response, open(outfile, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
            unzippath = args.type + "_HMMs"
            if not os.path.exists(unzippath):
                import zipfile
                zip_ref = zipfile.ZipFile(outfile, 'r')
                zip_ref.extractall(unzippath)

            if not os.path.exists(config["HMM_FOLDER"]):
                os.mkdir(config["HMM_FOLDER"])

            subdirs = os.listdir(unzippath)            
            for subd in subdirs:
                dlong = os.path.join(unzippath,subd)
                sublist = os.listdir(dlong)
                if 'HMM' in sublist:
                    for hmmfolder in os.listdir(os.path.join(dlong,'HMM')):
                        fromf=os.path.join(dlong,'HMM',hmmfolder)
                        print(fromf)
                        shutil.move(fromf,config["HMM_FOLDER"])
        print("%s HMMs installed. URL=%s Outfile=%s Dest=%s" % (args.type,url,outfile,config["HMM_FOLDER"]))
    if re.match("(hmm|search)",sys.argv[1].lower()):
        help = """
Usage:       PHYling %s <arguments>
version:     %s

Description: Script will search HMM set defined in config.txt against the genomes in pep
             
""" % ('search', version)
        arguments = sys.argv[2:]
        parser = argparse.ArgumentParser(description=help,add_help=True,formatter_class=argparse.RawDescriptionHelpFormatter)
        args = parser.parse_args(arguments)
        
        cmd = os.path.join(script_path, 'bin', 'phyling-01-hmmsearch.sh')
        logging.warning("running {}".format(cmd))
        subprocess.call(cmd)
