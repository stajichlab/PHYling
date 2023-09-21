from pathlib import Path
from io import StringIO
from Bio import Phylo

from phyling.phylotree import tree_generator


def get_newick(inputs):
    tree_generator_obj = tree_generator("upgma", 1, *inputs)
    tree = tree_generator_obj.get()
    tree_newick = StringIO()
    Phylo.write(tree, tree_newick, format="newick")
    tree_newick.seek(0)
    return tree_newick

def test_on_pep_single_tree():
    inputs = list(Path("tests/pep_align_single").iterdir())

    with open("tests/pep_single_upgma_tree.nw", "r") as f:
        assert get_newick(inputs).read() == f.read()


def test_on_pep_consensus_tree():
    inputs = list(Path("tests/pep_align_consensus").iterdir())

    with open("tests/pep_consensus_upgma_tree.nw", "r") as f:
        assert get_newick(inputs).read() == f.read()


def test_on_cds_single_tree():
    inputs = list(Path("tests/cds_align_single").iterdir())

    with open("tests/cds_single_upgma_tree.nw", "r") as f:
        assert get_newick(inputs).read() == f.read()


def test_on_cds_consensus_tree():
    inputs = list(Path("tests/cds_align_consensus").iterdir())

    with open("tests/cds_consensus_upgma_tree.nw", "r") as f:
        assert get_newick(inputs).read() == f.read()
