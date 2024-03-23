"""Some project specific configuration parameters."""

import time
from pathlib import Path


def runtime(start_time: float):
    """Measure the runtime of a process by given the start time."""
    elapsed_sec = time.monotonic() - start_time
    hours = int(elapsed_sec // 3600)
    mins = int((elapsed_sec // 60) % 60)
    secs = elapsed_sec % 60
    msg = f"{secs:.3f} seconds"
    msg = f"{mins} minutes " + msg if mins > 0 else msg
    msg = f"{hours} hours " + msg if hours > 0 else msg
    return msg


# download
database = "https://busco-data.ezlab.org/v5/data"
cfg_dir = Path.home() / ".phyling"
default_HMM = "HMM"  # default directory for HMM downloads

# libphyling
seqtype_dna = "DNA"
seqtype_pep = "peptide"
libphyling_checkpoint = ".align.ckp"
protein_ext = "faa"  # protein fasta files
cds_ext = "cds"  # coding sequence alignment file extension
aln_ext = "mfa"  # general alignment extension
prot_aln_ext = f"aa.{aln_ext}"  # protein alignment extension
cds_aln_ext = f"cds.{aln_ext}"  # cds alignment extension

# phylotree
phylotree_checkpoint = ".tree.ckp"
avail_tree_methods = {"ft": "FastTree", "raxml": "RAxML", "iqtree": "IQTree"}
selected_msas_dir = "selected_MSAs"
treeness_file = "treeness.tsv"
final_tree_file = "final_tree.nw"
concat_file_basename = "concat_alignments"
partition_ext = "partition"  # RAxML supported partition file extension
