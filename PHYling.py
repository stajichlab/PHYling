#!/usr/bin/env python3

#wrapper for PHYling

import sys, os, subprocess, inspect, tarfile, shutil, argparse
import urllib.request, configparser, re
import json, logging
from pathlib import Path

import zipfile

# TODO: general variables and paths need to be set and initialized
# here at the beginning

language = "en"

# user specific config file - the nane of this should be overrideable
Config_file = 'config.txt'

version = '0.1' #phyling release version

# sys.argv parsing should occur to get/set an alternative config file name

script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,script_path)

# this code loads up some config files which list where the HMM libraries
# can be downloaded from
URL_file      = os.path.join(script_path,"lib","urls.json")
Messages_file = os.path.join(script_path,"lib","messages.json") # abstract this later based on different languages
HMMs_URL = {}
Messages = {}
with open(URL_file,"r") as jsonfile:  
    jsonToPython = json.loads(jsonfile.read())
    HMMs_URL = jsonToPython['HMMs']

with open(Messages_file,"r") as jsonfile:  
    Messages = json.loads(jsonfile.read())[language]    


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
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # parse the rest of the argument
    args = parser.parse_args(arguments)
    cmd = os.path.join(script_path, 'bin', 'phyling-00-initialize.sh')
    subprocess.call(cmd)
elif subprog == 'download':
    help = Messages['commands']['download'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t','--type',default='fungi')
    args = parser.parse_args(arguments)
    url = ""
    outfile = args.type + "_HMMs.zip"

    #HEREHERE
    if args.type in HMMs_URL:
        url = HMMs_URL[args.type]

        if not os.path.exists(outfile):
            with urllib.request.urlopen(url) as response, open(outfile, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
                unzippath = args.type + "_HMMs"
            
            if not os.path.exists(unzippath):
                zip_ref = zipfile.ZipFile(outfile, 'r')
                zip_ref.extractall(unzippath)
            
            if not os.path.exists(config["HMM_FOLDER"]):
                os.mkdir(config["HMM_FOLDER"])

            for subdir in os.listdir(unzippath):
                dlong = os.path.join(unzippath,subdir)
                sublist = os.listdir(dlong)
                if 'HMM' in sublist:
                    for hmmfolder in os.listdir(os.path.join(dlong,'HMM')):
                        fromf=os.path.join(dlong,'HMM',hmmfolder)
                        print(fromf)
                        shutil.move(fromf,config["HMM_FOLDER"])
                        print("%s HMMs installed. URL=%s Outfile=%s Dest=%s" % (args.type,url,outfile,config["HMM_FOLDER"]))
        
elif re.match("(hmm|search)",subprog):
    help = Messages['commands']['search'] % ('search', version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args(arguments)
        
    cmd = os.path.join(script_path, 'bin', 'phyling-01-hmmsearch.sh')
    subprocess.call(cmd)

elif re.match("aln",subprog):        
    help = Messages['commands']['aln'] % ('aln', version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
    # hmmer or muscle for multiple alignment?
    args = parser.parse_args(arguments)
    # optional additional arguments might be -force to make sure to overwrite all files
    # and/or a 'clean' option to remove existing files first
    cmd = os.path.join(script_path, 'bin', 'phyling-02-makeunaln.sh')
    subprocess.call(cmd)
    cmd = os.path.join(script_path, 'bin', 'phyling-03-aln.sh')
    subprocess.call(cmd)

elif subprog == "concat":
    help = Messages['commands']['superaln'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args(arguments)

elif subprog == "phylo":
    help = Messages['commands']['phylo'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # hmmer or muscle for multiple alignment?
    args = parser.parse_args(arguments)

elif re.match(r"genetree",subprog):
    help = Messages['commands']['genetrees'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
    # hmmer or muscle for multiple alignment?
    args = parser.parse_args(arguments)
elif re.match(r"coal",subprog) or subprog == "astral":
    help = Messages['commands']['coalesce'] % (sys.argv[1], version)
    parser = argparse.ArgumentParser(description=help,add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args(arguments)

elif subprog == "version":
    print("PHYling v%s" % version)
elif subprog == "citation":
    print(Messages["citation"])
else:
    print("%s option not recognized" % sys.argv[1])
    print(default_help)
