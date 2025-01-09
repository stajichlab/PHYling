"""Phyling library."""

from enum import Enum

from .. import CFG_DIR


# General
class SeqTypes:
    DNA = "dna"
    PEP = "pep"
    RNA = "rna"


class FileExts:
    ALN = "mfa"  # general alignment extension
    PEP_ALN = f"aa.{ALN}"  # protein alignment extension
    CDS_ALN = f"cds.{ALN}"  # cds alignment extension


# Download
BUSCO_URL = "https://busco-data.ezlab.org/v5/data"
HMM_DIR = CFG_DIR / "HMM"  # default directory for HMM downloads
METADATA_FILE = CFG_DIR / "metadata.tsv"

# Align
ALIGN_METHODS = ("hmmalign", "muscle")


# Tree
class TreeMethods(Enum):
    FT = ("FastTree", ("fasttree", "Fasttree"))
    RAXML = ("RAxML-NG", ("raxml-ng",))
    IQTREE = ("IQTree", ("iqtree", "iqtree2"))

    def __init__(self, method: str, bins: tuple) -> None:
        self.method = method
        self.bins = bins


class TreeOutputFiles:
    MSAS_DIR = "selected_MSAs"
    TREENESS = "treeness.tsv"
    CONCAT = f"concat_alignments.{FileExts.ALN}"
    PARTITION = "concat_alignments.partition"
    TREE_NW = "final_tree.nw"
    TREE_IMG = "final_tree.png"
