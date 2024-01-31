"""Some project specific configuration parameters."""
from pathlib import Path

database = "https://busco-data.ezlab.org/v5/data"
cfg_dir = Path.home() / ".phyling"
default_HMM = "HMM"  # default directory for HMM downloads

prot_aln_ext = "aa.mfa"  # protein alignment extension
cds_aln_ext = "cds.mfa"  # cds alignment extension
protein_ext = "faa"  # protein fasta files
cds_ext = "cds"  # coding sequence alignment file extension
aln_ext = "mfa"  # general alignment extension

avail_tree_methods = {"upgma": "UPGMA", "nj": "Neighbor Joining", "ft": "FastTree"}
