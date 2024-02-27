from pathlib import Path

import pytest
from Bio import AlignIO, Phylo

from phyling.phylotree import concatenate_fasta, output_concat_file, mfa_to_finaltree

samples = [
    "Actinomucor_elegans_CBS_100.09",
    "Pilobolus_umbonatus_NRRL_6349",
    "Rhizopus_homothallicus_CBS_336.62",
    "Rhizopus_rhizopodiformis_NRRL_2570",
    "Zygorhynchus_heterogamous_NRRL_1489",
]


@pytest.fixture(scope="class")
def shared_tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("shared_tmpdir")


class TestFastaConcat:
    @pytest.fixture(scope="class", autouse=True)
    def prepare_data(self, shared_tmpdir):
        alignmentList = []
        self.__class__.inputs = list(Path("tests/pep_align").iterdir())
        for file in self.inputs:
            alignment = AlignIO.read(file, format="fasta")
            alignment.annotations["seqtype"] = "peptide"
            alignment.annotations["seqname"] = file.name
            alignmentList.append(alignment)
        self.__class__.concat_alignments = concatenate_fasta(samples, alignmentList, threads=1)
        self.__class__.concat_file, self.__class__.partition_file = output_concat_file(shared_tmpdir, self.concat_alignments)

    def test_concatenate_fasta_with_partition(self):
        print(self.concat_alignments)
        assert len(self.concat_alignments) == len(samples)
        assert len(self.concat_alignments.annotations["partition"]) == len(self.inputs)

    def test_output_concat_file(self):
        assert self.concat_file.is_file() is True and self.concat_file.stat().st_size > 0

    def test_output_partition_file(self):
        assert self.partition_file.is_file() is True and self.partition_file.stat().st_size > 0


class TestFastaConcatWoPartition(TestFastaConcat):
    @pytest.fixture(scope="class", autouse=True)
    def prepare_data(self, shared_tmpdir):
        alignmentList = []
        self.__class__.inputs = list(Path("tests/pep_align").iterdir())
        for file in self.inputs:
            alignment = AlignIO.read(file, format="fasta")
            alignmentList.append(alignment)
        self.__class__.concat_alignments = concatenate_fasta(samples, alignmentList, threads=1)
        self.__class__.concat_file, self.__class__.partition_file = output_concat_file(shared_tmpdir, self.concat_alignments)

    def test_concatenate_fasta_with_partition(self):
        print(self.concat_alignments)
        assert len(self.concat_alignments) == len(samples)

    def test_output_partition_file(self):
        assert self.partition_file is None


@pytest.mark.parametrize(
    "method, seqtype, concat",
    [
        ("upgma", "peptide", True),
        ("nj", "peptide", True),
        ("ft", "peptide", True),
        ("upgma", "peptide", False),
        ("nj", "peptide", False),
        ("ft", "peptide", False),
    ],
)
def test_mfa_to_finaltree_peptide(method, seqtype, concat, shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    final_tree = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method=method,
        seqtype=seqtype,
        samples=samples,
        concat=concat,
        top_n_toverr=2,
        threads=1,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True


@pytest.mark.parametrize(
    "method, seqtype, concat",
    [
        ("upgma", "DNA", True),
        ("nj", "DNA", True),
        ("ft", "DNA", True),
        ("upgma", "DNA", False),
        ("nj", "DNA", False),
        ("ft", "DNA", False),
    ],
)
def test_mfa_to_finaltree_cds(method, seqtype, concat, shared_tmpdir):
    inputs = list(Path("tests/cds_align").iterdir())
    final_tree = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method=method,
        seqtype=seqtype,
        samples=samples,
        concat=concat,
        top_n_toverr=2,
        threads=1,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True


def test_mfa_to_finaltree_top_n_toverr_out_of_range(shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    final_tree = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method="upgma",
        seqtype="peptide",
        samples=samples,
        concat=True,
        top_n_toverr=8,
        threads=1,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True


def test_mfa_to_finaltree_seqtype_KeyError(shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    with pytest.raises(KeyError):
        mfa_to_finaltree(
            inputs=inputs,
            output=shared_tmpdir,
            method="upgma",
            seqtype="RNA",
            samples=samples,
            concat=True,
            top_n_toverr=2,
            threads=1,
        )


def test_mfa_to_finaltree_method_KeyError(shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    with pytest.raises(KeyError):
        mfa_to_finaltree(
            inputs=inputs,
            output=shared_tmpdir,
            method="not_implemented_method",
            seqtype="peptide",
            samples=samples,
            concat=True,
            top_n_toverr=2,
            threads=1,
        )
