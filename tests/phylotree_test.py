from pathlib import Path

import pytest
from Bio import AlignIO, Phylo, SeqIO

import phyling.config
from phyling.phylotree import concatenate_fasta, tree_generator

samples = [
    "Actinomucor_elegans_CBS_100.09",
    "Pilobolus_umbonatus_NRRL_6349",
    "Rhizopus_homothallicus_CBS_336.62",
    "Rhizopus_rhizopodiformis_NRRL_2570",
    "Zygorhynchus_heterogamous_NRRL_1489",
]


@pytest.fixture(scope="module")
def shared_tmpdir(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp("shared_tmpdir")
    yield Path(tmpdir)


def concat(inputs, dir):
    alignmentList = []
    for file in inputs:
        alignmentList.append(AlignIO.read(file, format="fasta"))
    concat_alignments = concatenate_fasta(samples, alignmentList, threads=1)

    # Output the concatenated MSA to file
    output_concat = dir / f"concat_alignments.{phyling.config.aln_ext}"
    with open(output_concat, "w") as f:
        SeqIO.write(concat_alignments, f, format="fasta")
    inputs = [output_concat]
    return inputs


def test_on_pep_consensus():
    inputs = list(Path("tests/pep_align").iterdir())
    tree_generator_obj = tree_generator("peptide", "upgma", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_on_cds_consensus():
    inputs = list(Path("tests/cds_align").iterdir())
    tree_generator_obj = tree_generator("DNA", "upgma", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_on_pep_concat(tmp_path):
    inputs = list(Path("tests/pep_align").iterdir())
    inputs = concat(inputs, tmp_path)
    tree_generator_obj = tree_generator("peptide", "upgma", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_on_cds_concat(shared_tmpdir):
    inputs = list(Path("tests/cds_align").iterdir())
    inputs = concat(inputs, shared_tmpdir)
    tree_generator_obj = tree_generator("DNA", "upgma", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_nj_consensus():
    inputs = list(Path("tests/cds_align").iterdir())
    tree_generator_obj = tree_generator("DNA", "nj", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_nj_concat(shared_tmpdir):
    inputs = [shared_tmpdir / "concat_alignments.mfa"]
    tree_generator_obj = tree_generator("DNA", "nj", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_fasttree_consensus():
    inputs = list(Path("tests/cds_align").iterdir())
    tree_generator_obj = tree_generator("DNA", "ft", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True


def test_fasttree_concat(shared_tmpdir):
    inputs = [shared_tmpdir / "concat_alignments.mfa"]
    tree_generator_obj = tree_generator("DNA", "ft", *inputs)
    tree = tree_generator_obj.get(threads=1)
    assert isinstance(tree, Phylo.BaseTree.Tree) is True
