"""Some project specific configuration parameters."""
import os
from pathlib import Path

# CPUs
avail_cpus = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count()))

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
metadata = "metadata.tsv"

# libphyling
libphyling_checkpoint = ".align.ckp"
libphyling_precheck_params = {"inputs", "markerset", "markerset_cutoff", "evalue", "method", "non_trim"}
avail_align_methods = ("hmmalign", "muscle")

# phylotree
phylotree_checkpoint = ".tree.ckp"
phylotree_precheck_params = {"method", "top_n_toverr", "concat", "partition", "inputs"}
avail_tree_methods = {"ft": "FastTree", "raxml": "RAxML", "iqtree": "IQTree"}
