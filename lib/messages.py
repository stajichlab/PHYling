#!/usr/bin/env python3
import json,os,inspect

Messages = { 
    "en": {
        "commands": {
            "default": """
Usage:       PHYling <command> <arguments>
version:     %s
Description: PHYling is a package of scripts to extract phylogenomic markers
Dependencies:  HMMER3 MUSCLE RAxML

init        Setup or recreate initialize files
search      Run HMMsearch of markers against proteomes
aln         Construct individual gene alignments and trim
superaln    Concatenate gene alignments into superalignment
phylo       Run phylogenetic analysis on superalignment
genetrees   Run phylogenetic analysis on individual genes
coalesce    Run ASTRAL coalescence of gene trees

download    Download HMM markers


Written by Jason Stajich (2014-2017)
jason.stajich[at]ucr.edu or jasonstajich.phd[AT]gmail.com

Initially written https://github.com/1KFG/Phylogenomics and
https://github.com/stajichlab/phyling
""",
            "initialize": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Script will setup initial data directory and files. Expects an HMM folder
             config.txt and 'pep' folder as specified in config.txt file. An optional 'cds' folder
""",
        "download": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Download HMM files for pre-defined phylogenomic markers. 

Arguments:   -t fungi [default]
""",
            "search": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Search HMM set defined in config.txt against the genomes in pep
             
""",
            "aln": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Construct unaligned fasta files of protein and cds (if found) and
             perform multiple alignments
             
""",
            "superaln": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Concatenate gene alignments into a superalignment

Options:    -a [denovo|hmmalign]
             
""",        
            "phylo": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Run Phylogenetic reconstruction on superalignment

Options:     -t [raxml|iqtree|fasttree]
             -m [aa|cds]
             -a [denovo|hmmalign]
""",
            "genetrees": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Run Phylogenetic reconstruction on each gene from marker

Options:     -t [raxml|fasttree]
             -m [aa|cds]
""",
            "coalesce": """
Usage:       PHYling %s <arguments>
version:     %s

Description: Run ASTRAL gene tree reconciliation software

Options:     -t [raxml|iqtree|fasttree] - source tree names
             -m [aa|cds]                - protein or cds trees
             -b                         - sample bootstrap files if generated (RAXML)
"""

        },
        "citation": """
Stajich JE. 2017 PHYling: Phylogenomic pipeline from core markers on genomes to unassembled reads. https://github.com/stajichlab/PHYling_unified
"""
    }
}
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

outfile = os.path.join(script_path,"messages.json")
#print(outfile)
of = open(outfile, "w")
of.write(json.dumps(Messages))
