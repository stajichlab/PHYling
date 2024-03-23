from pathlib import Path

import pytest
from Bio import Phylo

import phyling.config
from phyling.phylotree import mfa_to_finaltree, phylotree

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


@pytest.mark.parametrize(
    "method, concat, partition",
    [
        ("ft", True, None),
        ("raxml", True, None),
        ("iqtree", True, None),
        ("ft", True, "seq"),
        ("raxml", True, "seq"),
        ("iqtree", True, "seq"),
        ("ft", False, None),
        ("raxml", False, None),
        ("iqtree", False, None),
    ],
)
def test_mfa_to_finaltree_peptide(method, concat, partition, shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    final_tree, _ = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method=method,
        seqtype=phyling.config.seqtype_pep,
        samples=samples,
        concat=concat,
        top_n_toverr=5,
        partition=partition,
        threads=1,
        rerun=0,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True
    if concat:
        assert (Path(shared_tmpdir) / f"{phyling.config.concat_file_basename}.{phyling.config.aln_ext}").is_file()
        if partition:
            assert (Path(shared_tmpdir) / f"{phyling.config.concat_file_basename}.{phyling.config.partition_ext}").is_file()


@pytest.mark.parametrize(
    "method, concat, partition",
    [
        ("ft", True, None),
        ("raxml", True, None),
        ("iqtree", True, None),
        ("ft", True, "seq"),
        ("raxml", True, "seq"),
        ("iqtree", True, "seq"),
        ("ft", True, "codon"),
        ("raxml", True, "codon"),
        ("iqtree", True, "codon"),
        ("ft", True, "seq+codon"),
        ("raxml", True, "seq+codon"),
        ("iqtree", True, "seq+codon"),
        ("ft", False, None),
        ("raxml", False, None),
        ("iqtree", False, None),
    ],
)
def test_mfa_to_finaltree_cds(method, concat, partition, shared_tmpdir):
    inputs = list(Path("tests/cds_align").iterdir())
    final_tree, _ = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method=method,
        seqtype=phyling.config.seqtype_dna,
        samples=samples,
        concat=concat,
        top_n_toverr=5,
        partition=partition,
        threads=1,
        rerun=0,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True
    if concat:
        assert (Path(shared_tmpdir) / f"{phyling.config.concat_file_basename}.{phyling.config.aln_ext}").is_file()
        if partition:
            assert (Path(shared_tmpdir) / f"{phyling.config.concat_file_basename}.{phyling.config.partition_ext}").is_file()


@pytest.mark.parametrize(
    "inputs, seqtype",
    [
        (list(Path("tests/pep_align").iterdir()), "peptide"),
        (list(Path("tests/cds_align").iterdir()), "DNA"),
    ],
)
def test_mfa_to_finaltree_top_n_toverr_out_of_range(inputs, seqtype, shared_tmpdir):
    final_tree, tree_gen_obj = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method="ft",
        seqtype=seqtype,
        samples=samples,
        concat=True,
        top_n_toverr=8,
        partition=None,
        threads=1,
        rerun=0,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True
    assert isinstance(tree_gen_obj.get_mfa(), list) is True
    assert len(tree_gen_obj.get_mfa()) == 5


@pytest.mark.parametrize(
    "inputs, seqtype",
    [
        ([Path("tests/pep_align/112646at4751.aa.mfa")], "peptide"),
        ([Path("tests/cds_align/112646at4751.cds.mfa")], "DNA"),
    ],
)
def test_mfa_to_finaltree_single_input(inputs, seqtype, shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    final_tree, tree_gen_obj = mfa_to_finaltree(
        inputs=inputs,
        output=shared_tmpdir,
        method="ft",
        seqtype="peptide",
        samples=samples,
        concat=True,
        top_n_toverr=5,
        partition=None,
        threads=1,
        rerun=0,
    )
    assert isinstance(final_tree, Phylo.BaseTree.Tree) is True


def test_mfa_to_finaltree_seqtype_KeyError(shared_tmpdir):
    inputs = list(Path("tests/pep_align").iterdir())
    with pytest.raises(KeyError):
        mfa_to_finaltree(
            inputs=inputs,
            output=shared_tmpdir,
            method="ft",
            seqtype="RNA",
            samples=samples,
            concat=True,
            top_n_toverr=2,
            partition=None,
            threads=1,
            rerun=0,
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
            partition=None,
            threads=1,
            rerun=0,
        )


class TestCheckpointRerun:
    @pytest.fixture(scope="class", autouse=True)
    def filepath(self, shared_tmpdir):
        self.__class__.checkpoint = shared_tmpdir / phyling.config.phylotree_checkpoint
        self.__class__.final_tree_file = shared_tmpdir / phyling.config.final_tree_file
        self.__class__.treeness_file = shared_tmpdir / phyling.config.treeness_file
        self.__class__.selected_msas_dir = shared_tmpdir / phyling.config.selected_msas_dir
        self.__class__.concat_file = shared_tmpdir / f"{phyling.config.concat_file_basename}.{phyling.config.aln_ext}"
        self.__class__.partition_file = shared_tmpdir / f"{phyling.config.concat_file_basename}.{phyling.config.partition_ext}"

    def _run_phylotree(self, output, method, top_n_toverr, concat, partition):
        phylotree(
            inputs=None,
            input_dir=Path("tests/pep_align"),
            output=output,
            method=method,
            top_n_toverr=top_n_toverr,
            concat=concat,
            partition=partition,
            figure=False,
            threads=1,
        )

    def _get_mtime(self):
        checkpoint = self.checkpoint.stat().st_mtime_ns
        final_tree_file = self.final_tree_file.stat().st_mtime_ns
        concat_file = self.concat_file.stat().st_mtime_ns if self.concat_file.exists() else False
        partition_file = self.partition_file.stat().st_mtime_ns if self.partition_file.exists() else False
        treeness_file = self.treeness_file.stat().st_mtime_ns
        selected_msas_dir = self.selected_msas_dir.stat().st_mtime_ns
        return checkpoint, final_tree_file, concat_file, partition_file, treeness_file, selected_msas_dir

    def test_checkpoint_creation(self, shared_tmpdir):
        assert not self.checkpoint.exists()
        self._run_phylotree(output=shared_tmpdir, method="ft", top_n_toverr=5, concat=True, partition="seq")
        assert self.checkpoint.is_file()

    def test_no_rerun(self, shared_tmpdir):
        checkpoint_before = self.checkpoint.stat().st_mtime_ns
        with pytest.raises(SystemExit):
            self._run_phylotree(output=shared_tmpdir, method="ft", top_n_toverr=5, concat=True, partition="seq")
        assert checkpoint_before == self.checkpoint.stat().st_mtime_ns

    def test_input_changed(self, shared_tmpdir):
        checkpoint_before = self.checkpoint.stat().st_mtime_ns
        with pytest.raises(SystemExit):
            phylotree(
                inputs=None,
                input_dir=Path("tests/cds_align"),
                output=shared_tmpdir,
                method="ft",
                top_n_toverr=5,
                concat=True,
                partition="seq",
                figure=False,
                threads=1,
            )
        assert checkpoint_before == self.checkpoint.stat().st_mtime_ns

    @pytest.mark.parametrize(
        "method, top_n_toverr, concat, partition",
        [
            ("raxml", 5, True, "seq"),
            ("iqtree", 5, True, "seq"),
            ("ft", 5, True, "seq"),
        ],
    )
    def test_rerun_final_tree(self, method, top_n_toverr, concat, partition, shared_tmpdir):
        (
            checkpoint_before,
            final_tree_file_before,
            concat_file_before,
            partition_file_before,
            treeness_file_before,
            selected_msas_dir_before,
        ) = self._get_mtime()

        self._run_phylotree(output=shared_tmpdir, method=method, top_n_toverr=top_n_toverr, concat=concat, partition=partition)

        assert checkpoint_before != self.checkpoint.stat().st_mtime_ns
        assert final_tree_file_before != self.final_tree_file.stat().st_mtime_ns
        assert concat_file_before == self.concat_file.stat().st_mtime_ns
        assert partition_file_before == self.partition_file.stat().st_mtime_ns
        assert treeness_file_before == self.treeness_file.stat().st_mtime_ns
        assert selected_msas_dir_before == self.selected_msas_dir.stat().st_mtime_ns

    @pytest.mark.parametrize(
        "method, top_n_toverr, concat, partition",
        [
            ("ft", 5, True, None),
            ("raxml", 5, True, "seq"),
        ],
    )
    def test_rerun_concat(self, method, top_n_toverr, concat, partition, shared_tmpdir):
        (
            checkpoint_before,
            final_tree_file_before,
            concat_file_before,
            partition_file_before,
            treeness_file_before,
            selected_msas_dir_before,
        ) = self._get_mtime()

        self._run_phylotree(output=shared_tmpdir, method=method, top_n_toverr=top_n_toverr, concat=concat, partition=partition)

        assert checkpoint_before != self.checkpoint.stat().st_mtime_ns
        assert final_tree_file_before != self.final_tree_file.stat().st_mtime_ns
        assert concat_file_before != self.concat_file.stat().st_mtime_ns
        assert partition_file_before != (self.partition_file.stat().st_mtime_ns if self.partition_file.exists() else False)
        assert treeness_file_before == self.treeness_file.stat().st_mtime_ns
        assert selected_msas_dir_before == self.selected_msas_dir.stat().st_mtime_ns

    @pytest.mark.parametrize(
        "method, top_n_toverr, concat, partition",
        [
            ("raxml", 4, True, "seq"),
            ("ft", 5, True, None),
        ],
    )
    def test_rerun_msa_selection(self, method, top_n_toverr, concat, partition, shared_tmpdir):
        (
            checkpoint_before,
            final_tree_file_before,
            concat_file_before,
            partition_file_before,
            treeness_file_before,
            selected_msas_dir_before,
        ) = self._get_mtime()

        self._run_phylotree(output=shared_tmpdir, method=method, top_n_toverr=top_n_toverr, concat=concat, partition=partition)

        assert checkpoint_before != self.checkpoint.stat().st_mtime_ns
        assert final_tree_file_before != self.final_tree_file.stat().st_mtime_ns
        assert concat_file_before != self.concat_file.stat().st_mtime_ns
        assert partition_file_before != (self.partition_file.stat().st_mtime_ns if self.partition_file.exists() else False)
        assert treeness_file_before != self.treeness_file.stat().st_mtime_ns
        assert selected_msas_dir_before != self.selected_msas_dir.stat().st_mtime_ns

    @pytest.mark.parametrize(
        "method, top_n_toverr, concat, partition",
        [
            ("ft", 4, False, None),
            ("ft", 5, True, None),
            ("ft", 5, False, None),
        ],
    )
    def test_switch_type(self, method, top_n_toverr, concat, partition, shared_tmpdir):
        (
            checkpoint_before,
            final_tree_file_before,
            concat_file_before,
            partition_file_before,
            treeness_file_before,
            selected_msas_dir_before,
        ) = self._get_mtime()

        self._run_phylotree(output=shared_tmpdir, method=method, top_n_toverr=top_n_toverr, concat=concat, partition=partition)

        assert checkpoint_before != self.checkpoint.stat().st_mtime_ns
        assert final_tree_file_before != self.final_tree_file.stat().st_mtime_ns
        assert concat_file_before != (self.concat_file.stat().st_mtime_ns if self.concat_file.exists() else False)
        assert partition_file_before == self.partition_file.exists() is False
        assert treeness_file_before != self.treeness_file.stat().st_mtime_ns
        assert selected_msas_dir_before != self.selected_msas_dir.stat().st_mtime_ns

    @pytest.mark.parametrize(
        "method, top_n_toverr, concat, partition",
        [
            ("raxml", 5, False, None),
            ("iqtree", 4, False, None),
            ("ft", 5, False, None),
        ],
    )
    def test_consensus_rerun(self, method, top_n_toverr, concat, partition, shared_tmpdir):
        (
            checkpoint_before,
            final_tree_file_before,
            concat_file_before,
            partition_file_before,
            treeness_file_before,
            selected_msas_dir_before,
        ) = self._get_mtime()

        self._run_phylotree(output=shared_tmpdir, method=method, top_n_toverr=top_n_toverr, concat=concat, partition=partition)

        assert checkpoint_before != self.checkpoint.stat().st_mtime_ns
        assert final_tree_file_before != self.final_tree_file.stat().st_mtime_ns
        assert concat_file_before == self.concat_file.exists() is False
        assert partition_file_before == self.partition_file.exists() is False
        assert treeness_file_before != self.treeness_file.stat().st_mtime_ns
        assert selected_msas_dir_before != self.selected_msas_dir.stat().st_mtime_ns
