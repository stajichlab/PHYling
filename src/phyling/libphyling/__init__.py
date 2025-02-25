"""Phyling library."""

from enum import Enum


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
METADATA_FILE = ".metadata"

# Align
ALIGN_METHODS = ("hmmalign", "muscle")


# Tree
class TreeMethods(Enum):
    FT = ("FastTree", ("VeryFastTree",), ("JC", "GTR"), ("JTT", "WAG", "LG"))
    RAXML = ("RAxML-NG", ("raxml-ng",), ("raxml",), ("raxml",))
    IQTREE = ("IQTree", ("iqtree", "iqtree2"), (), ())

    def __init__(self, method: str, bins: tuple, dna_model: tuple, pep_model: tuple) -> None:
        self.method = method
        self.bins = bins
        self.dna_model = dna_model
        self.pep_model = pep_model


class TreeOutputFiles:
    MSAS_DIR = "selected_MSAs"
    TREENESS = "treeness.tsv"
    CONCAT = f"concat_alignments.{FileExts.ALN}"
    PARTITION = "concat_alignments.partition"
    TREE_NW = "final_tree.nw"
    TREE_IMG = "final_tree.png"
