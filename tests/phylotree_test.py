from pathlib import Path
from Bio import Phylo

from phyling.phylotree import tree_generator


def test_on_pep_single_tree():
    inputs = list(Path("tests/pep_align_single").iterdir())
    tree_generator_obj = tree_generator("upgma", 1, *inputs)
    tree = tree_generator_obj.get()
    assert isinstance(tree, Phylo.BaseTree.Tree) == True


# def test_on_pep_consensus_tree():
#     inputs = list(Path("tests/pep_align_consensus").iterdir())
#     tree_generator_obj = tree_generator("upgma", 1, *inputs)
#     tree = tree_generator_obj.get()
#     assert isinstance(tree, Phylo.BaseTree.Tree) == True


def test_on_cds_single_tree():
    inputs = list(Path("tests/cds_align_single").iterdir())
    tree_generator_obj = tree_generator("upgma", 1, *inputs)
    tree = tree_generator_obj.get()
    assert isinstance(tree, Phylo.BaseTree.Tree) == True


# def test_on_cds_consensus_tree():
#     inputs = list(Path("tests/cds_align_consensus").iterdir())
#     tree_generator_obj = tree_generator("upgma", 1, *inputs)
#     tree = tree_generator_obj.get()
#     assert isinstance(tree, Phylo.BaseTree.Tree) == True


# def test_nj():
#     inputs = list(Path("tests/cds_align_consensus").iterdir())
#     tree_generator_obj = tree_generator("nj", 1, *inputs)
#     tree = tree_generator_obj.get()
#     assert isinstance(tree, Phylo.BaseTree.Tree) == True


# def test_fasttree():
#     inputs = list(Path("tests/cds_align_consensus").iterdir())
#     tree_generator_obj = tree_generator("ft", 1, *inputs)
#     tree = tree_generator_obj.get()
#     assert isinstance(tree, Phylo.BaseTree.Tree) == True
