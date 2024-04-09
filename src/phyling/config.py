"""Some project specific configuration parameters."""

from pathlib import Path

# Common sequence and file constants
seqtype_cds = "cds"
seqtype_pep = "pep"
protein_ext = "faa"  # protein fasta files
cds_ext = "cds"  # coding sequence alignment file extension
aln_ext = "mfa"  # general alignment extension
prot_aln_ext = f"aa.{aln_ext}"  # protein alignment extension
cds_aln_ext = f"cds.{aln_ext}"  # cds alignment extension

# download
database = "https://busco-data.ezlab.org/v5/data"
cfg_dir = Path.home() / ".phyling"
default_HMM = "HMM"  # default directory for HMM downloads

# libphyling
libphyling_checkpoint = ".align.ckp"

# phylotree
phylotree_checkpoint = ".tree.ckp"
avail_tree_methods = {"ft": "FastTree", "raxml": "RAxML", "iqtree": "IQTree"}
