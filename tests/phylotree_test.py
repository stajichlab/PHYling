from __future__ import annotations

import logging
import shutil
from copy import deepcopy
from itertools import permutations, product
from pathlib import Path

import pytest
from Bio.AlignIO import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree

import phyling._internal._config as _config
import phyling.exception as exception
from phyling._internal._libtree import MFA2Tree, MFA2TreeWrapper, OutputPrecheck, determine_samples_and_seqtype


class TestMFA2Tree:
    @pytest.mark.parametrize(
        "file, partition, seqtype",
        (
            ("tests/data/pep_align/10at10240.aa.mfa", None, "pep"),
            ("tests/data/cds_align/14at10240.cds.mfa", None, "cds"),
            ("tests/data/pep_align/10at10240.aa.mfa", None, None),
            ("tests/data/cds_align/14at10240.cds.mfa", None, None),
            ("tests/data/concat_pep.mfa", "tests/data/concat_pep.partition", "pep"),
            ("tests/data/concat_cds.mfa", "tests/data/concat_cds.partition", "cds"),
        ),
    )
    def test_init(self, file: str, partition: str, seqtype: str):
        m = MFA2Tree(file, partition, seqtype=seqtype)
        assert m.name == Path(file).name
        assert isinstance(m.msa, MultipleSeqAlignment)

    def test_eq_lt(self):
        a = MFA2Tree("tests/data/pep_align/10at10240.aa.mfa", seqtype="pep")  # 1.3572
        b = MFA2Tree("tests/data/pep_align/10at10240.aa.mfa", seqtype="pep")
        c = MFA2Tree("tests/data/pep_align/74at10240.aa.mfa", seqtype="pep")  # 1.5419
        for x in (a, b, c):
            x.build("ft")
            x.compute_toverr()
        assert a == b < c

    def test_eq_lt_attributeerror(self):
        a = MFA2Tree("tests/data/pep_align/10at10240.aa.mfa", seqtype="pep")  # 1.3572
        b = MFA2Tree("tests/data/pep_align/10at10240.aa.mfa", seqtype="pep")
        c = MFA2Tree("tests/data/pep_align/74at10240.aa.mfa", seqtype="pep")  # 1.5419
        with pytest.raises(AttributeError, match=r"build\(\) need to be run first"):
            a == b
        with pytest.raises(AttributeError, match=r"build\(\) need to be run first"):
            a < c
        for x in (a, b, c):
            x.build("ft")
        with pytest.raises(AttributeError, match=r"compute_toverr\(\) need to be run first"):
            a == b
        with pytest.raises(AttributeError, match=r"compute_toverr\(\) need to be run first"):
            a < c

    def test_eq_lt_seqtypeerror(self):
        a = MFA2Tree("tests/data/pep_align/10at10240.aa.mfa", seqtype="pep")  # 1.3572
        b = MFA2Tree("tests/data/cds_align/10at10240.cds.mfa", seqtype="cds")  # 2.1163
        c = MFA2Tree("tests/data/pep_align/74at10240.aa.mfa", seqtype="pep")  # 1.3572
        for x in (a, b, c):
            x.build("ft")
            x.compute_toverr()
        with pytest.raises(exception.SeqtypeError, match="Items represent different seqtypes"):
            a < b
        with pytest.raises(exception.SeqtypeError, match="Items represent different seqtypes"):
            b > c

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "input, partition, seqtype, method",
        [
            (input, partition, seqtype, method)
            for (input, partition, seqtype), method in product(
                (
                    ("tests/data/pep_align/10at10240.aa.mfa", None, "pep"),
                    ("tests/data/cds_align/10at10240.cds.mfa", None, "cds"),
                    ("tests/data/concat_pep.mfa", "tests/data/concat_pep.partition", "pep"),
                    ("tests/data/concat_cds.mfa", "tests/data/concat_cds.partition", "cds"),
                ),
                ("ft", "raxml", "iqtree"),
            )
        ],
    )
    def test_build(self, input, partition, seqtype, method):
        m = MFA2Tree(input, partition, seqtype=seqtype)
        m.build(method)
        assert m.method == method
        assert isinstance(m.tree, Tree)

    @pytest.mark.parametrize(
        "file, seqtype",
        (
            ("tests/data/pep_align/10at10240.aa.mfa", "pep"),
            ("tests/data/cds_align/14at10240.cds.mfa", "cds"),
        ),
    )
    def test_compute_toverr(self, file: str, seqtype: str):
        m = MFA2Tree(file, seqtype=seqtype)  # 1.3572
        m.build("ft")
        m.compute_toverr()
        assert isinstance(m.toverr, float)


class TestMFA2TreeWrapper:
    files_pep = (
        "tests/data/pep_align/10at10240.aa.mfa",
        "tests/data/pep_align/14at10240.aa.mfa",
        "tests/data/pep_align/22at10240.aa.mfa",
    )

    files_cds = (
        "tests/data/cds_align/14at10240.cds.mfa",
        "tests/data/cds_align/10at10240.cds.mfa",
        "tests/data/cds_align/22at10240.cds.mfa",
    )

    wrapper_pep = MFA2TreeWrapper(files_pep, seqtype="pep")
    wrapper_cds = MFA2TreeWrapper(files_cds, seqtype="cds")
    built_pep, built_cds = deepcopy(wrapper_pep), deepcopy(wrapper_cds)
    for x in (built_pep, built_cds):
        x.build("ft")
    sorted_pep, sorted_cds = deepcopy(built_pep), deepcopy(built_cds)
    for x in (sorted_pep, sorted_cds):
        x.compute_toverr()

    @pytest.mark.parametrize("wrapper, seqtype", ((built_pep, "pep"), (built_cds, "cds")))
    def test_init(self, wrapper: MFA2TreeWrapper, seqtype: str):
        assert wrapper.seqtype == seqtype
        assert wrapper.is_sorted is False

    @pytest.mark.parametrize("wrapper, sorted_wrapper", ((wrapper_pep, sorted_pep), (wrapper_cds, sorted_cds)))
    def test_getitem(self, wrapper: MFA2TreeWrapper, sorted_wrapper: MFA2TreeWrapper):
        assert len(wrapper) == 3
        assert all(isinstance(x, MFA2Tree) for x in wrapper)
        w_sub = sorted_wrapper[0:2]
        assert isinstance(w_sub, MFA2TreeWrapper)
        assert len(w_sub) == 2
        assert all(isinstance(x, MFA2Tree) for x in w_sub)

    @pytest.mark.parametrize("wrapper, built_wrapper", ((wrapper_pep, built_pep), (wrapper_cds, built_cds)))
    def test_getitem_attributeerror(self, wrapper: MFA2TreeWrapper, built_wrapper: MFA2TreeWrapper):
        with pytest.raises(AttributeError, match=r"build\(\) need to be run first"):
            wrapper[0:2]
        with pytest.raises(AttributeError, match=r"compute_toverr\(\) need to be run first"):
            built_wrapper[0:2]

    @pytest.mark.parametrize("wrapper", (built_pep, built_cds))
    def test_sort(self, wrapper: MFA2TreeWrapper):
        w = deepcopy(wrapper)
        w.compute_toverr()
        w.sort()
        assert w.is_sorted is True
        assert w[0] > w[1] > w[2]

    @pytest.mark.parametrize("wrapper, built_wrapper", ((wrapper_pep, built_pep), (wrapper_cds, built_cds)))
    def test_sort_attributeerror(self, wrapper: MFA2TreeWrapper, built_wrapper: MFA2TreeWrapper):
        with pytest.raises(AttributeError, match=r"build\(\) need to be run first"):
            wrapper.sort()
        with pytest.raises(AttributeError, match=r"compute_toverr\(\) need to be run first"):
            built_wrapper.sort()

    @pytest.mark.slow
    @pytest.mark.parametrize("wrapper, method", list(product((wrapper_pep, wrapper_cds), ("ft", "raxml", "iqtree"))))
    def test_build(self, wrapper: MFA2TreeWrapper, method: str):
        w = deepcopy(wrapper)
        w.build(method)
        assert w.method == method

    @pytest.mark.parametrize("wrapper", (built_pep, built_cds))
    def test_compute_toverr(self, wrapper: MFA2TreeWrapper):
        w = deepcopy(wrapper)
        w.compute_toverr()
        assert all(isinstance(x.toverr, float) for x in w)

    @pytest.mark.parametrize("wrapper", (sorted_pep, sorted_cds))
    def test_get_consensus_tree(self, wrapper: MFA2TreeWrapper):
        w = deepcopy(wrapper)
        assert isinstance(w.get_consensus_tree(), Tree)

    @pytest.mark.parametrize(
        "wrapper, samples, partition",
        (
            (sorted_pep, None, None),
            (sorted_cds, None, None),
            (sorted_pep, None, "seq"),
            (sorted_cds, None, "seq"),
            (sorted_cds, None, "codon"),
            (sorted_cds, None, "seq+codon"),
            (
                sorted_pep,
                ("Anomala_cuprea_entomopoxvirus", "Canarypox_virus", "Cowpox_virus", "Goatpox_virus", "Variola_virus"),
                "seq",
            ),
        ),
    )
    def test_concat(self, tmp_path: Path, wrapper: MFA2TreeWrapper, samples: tuple, partition: str):
        w = deepcopy(wrapper)
        r = w.concat(tmp_path, samples, partition=partition)
        assert (tmp_path / f"concat_alignments.{_config.aln_ext}").is_file() is True
        assert isinstance(r, MFA2Tree)
        assert r.name == f"concat_alignments.{_config.aln_ext}"
        if partition:
            assert (tmp_path / "concat_alignments.partition").is_file() is True


class TestOutputPrecheck:
    @pytest.fixture(scope="function")
    def copy_phylotree_ckp(self, tmp_path: Path) -> Path:
        shutil.copy("tests/data/tree/.tree.ckp", tmp_path / ".tree.ckp")
        OutputPrecheck.setup(folder=tmp_path)
        return tmp_path

    def test_precheck_new(self, tmp_path: Path):
        assert (tmp_path / "output").exists() is False
        OutputPrecheck.setup(folder=tmp_path / "output")
        params = {
            "method": "ft",
            "top_n_toverr": 50,
            "concat": False,
            "partition": None,
            "inputs": "fakechecksum",
        }
        r, _ = OutputPrecheck.precheck(params)
        assert r == 0
        assert (tmp_path / "output").is_dir() is True

    def test_precheck_folder_exists(self, tmp_path: Path):
        OutputPrecheck.setup(folder=tmp_path)
        params = {
            "method": "ft",
            "top_n_toverr": 18,
            "concat": False,
            "partition": None,
            "inputs": "fakechecksum",
        }
        r, _ = OutputPrecheck.precheck(params)
        assert r == 0
        assert tmp_path.is_dir() is True

    def test_precheck_notadirectoryerror(self, tmp_path: Path):
        (tmp_path / "output").touch()
        OutputPrecheck.setup(folder=tmp_path / "output")
        params = {
            "method": "ft",
            "top_n_toverr": 18,
            "concat": False,
            "partition": None,
            "inputs": "fakechecksum",
        }
        with pytest.raises(NotADirectoryError, match="already existed but not a folder"):
            OutputPrecheck.precheck(params)

    @pytest.mark.parametrize("method1, method2", list(permutations(("ft", "raxml", "iqtree"), 2)))
    def test_precheck_rerun_status_3(self, method1: str, method2: str, tmp_path: Path):
        OutputPrecheck.setup(folder=tmp_path)
        msas_dir = tmp_path / "selected_MSAs"
        msas_dir.mkdir()
        for file in Path("tests/data/pep_align").glob("*.mfa"):
            (msas_dir / file.name).symlink_to(file.absolute())

        prev_params = {
            "method": method1,
            "top_n_toverr": 18,
            "concat": True,
            "partition": None,
            "inputs": "fakechecksum",
        }
        fake_input = tmp_path / "input"
        fake_input.touch()
        OutputPrecheck.save_checkpoint(prev_params, MFA2TreeWrapper([fake_input], seqtype="pep"))
        cur_params = {
            "method": method2,
            "top_n_toverr": 18,
            "concat": True,
            "partition": None,
            "inputs": "fakechecksum",
        }
        r, _ = OutputPrecheck.precheck((cur_params))
        assert r == 3

    @pytest.mark.parametrize(
        "method, partition",
        tuple(
            product(
                tuple(product(("ft", "raxml", "iqtree"), repeat=2)),
                tuple(permutations((None, "seq", "codon", "seq+codon"), 2)),
            )
        ),
    )
    def test_precheck_rerun_status_2(self, method: str, partition: str, tmp_path: Path):
        OutputPrecheck.setup(folder=tmp_path)
        msas_dir = tmp_path / "selected_MSAs"
        msas_dir.mkdir()
        for file in Path("tests/data/pep_align").glob("*.mfa"):
            (msas_dir / file.name).symlink_to(file.absolute())

        prev_params = {
            "method": method[0],
            "top_n_toverr": 18,
            "concat": True,
            "partition": partition[0],
            "inputs": "fakechecksum",
        }
        fake_input = tmp_path / "input"
        fake_input.touch()
        OutputPrecheck.save_checkpoint(prev_params, MFA2TreeWrapper([fake_input], seqtype="pep"))
        cur_params = {
            "method": method[1],
            "top_n_toverr": 18,
            "concat": True,
            "partition": partition[1],
            "inputs": "fakechecksum",
        }
        r, _ = OutputPrecheck.precheck(cur_params)
        assert r == 2

    @pytest.mark.parametrize(
        "top_n, partition",
        tuple(
            product(
                ((18, 5),),
                tuple(product(("ft", "raxml", "iqtree"), repeat=2)),
            )
        ),
    )
    def test_precheck_rerun_status_1(self, top_n: int, partition: str, tmp_path: Path):
        OutputPrecheck.setup(folder=tmp_path)
        msas_dir = tmp_path / "selected_MSAs"
        msas_dir.mkdir()
        for file in Path("tests/data/pep_align").glob("*.mfa"):
            (msas_dir / file.name).symlink_to(file.absolute())

        prev_params = {
            "method": "ft",
            "top_n_toverr": top_n[0],
            "concat": True,
            "partition": partition[0],
            "inputs": "fakechecksum",
        }
        fake_input = tmp_path / "input"
        fake_input.touch()
        OutputPrecheck.save_checkpoint(prev_params, MFA2TreeWrapper([fake_input], seqtype="pep"))
        cur_params = {
            "method": "ft",
            "top_n_toverr": top_n[1],
            "concat": True,
            "partition": partition[1],
            "inputs": "fakechecksum",
        }
        r, _ = OutputPrecheck.precheck(cur_params)
        assert r == 1

    @pytest.mark.parametrize(
        "method, top_n, concat, partition",
        list(
            product(
                (("ft", "ft"), ("ft", "raxml")),
                ((10, 10), (10, 50)),
                ((False, False), (True, False), (False, True)),
                tuple(product((None, "seq", "codon", "seq+codon"), repeat=2)),
            )
        ),
    )
    def test_precheck_rerun_status_0(self, method: str, top_n: int, concat: bool, partition: str, tmp_path: Path):
        OutputPrecheck.setup(folder=tmp_path)
        msas_dir = tmp_path / "selected_MSAs"
        msas_dir.mkdir()
        for file in Path("tests/data/pep_align").glob("*.mfa"):
            (msas_dir / file.name).symlink_to(file.absolute())

        prev_params = {
            "method": method[0],
            "top_n_toverr": top_n[0],
            "concat": concat[0],
            "partition": partition[0],
            "inputs": "fakechecksum",
        }
        fake_input = tmp_path / "input"
        fake_input.touch()
        OutputPrecheck.save_checkpoint(prev_params.copy(), MFA2TreeWrapper([fake_input], seqtype="pep"))
        cur_params = {
            "method": method[1],
            "top_n_toverr": top_n[1],
            "concat": concat[1],
            "partition": partition[1],
            "inputs": "fakechecksum",
        }
        print(prev_params)
        if prev_params == cur_params:
            with pytest.raises(SystemExit, match="Files not changed and parameters are identical to the previous run"):
                OutputPrecheck.precheck(cur_params)
        else:
            r, _ = OutputPrecheck.precheck(cur_params)
            assert r == 0

    def test_load_checkpoint(self):
        OutputPrecheck.setup(folder=Path("tests/data/tree"))
        params, wrapper = OutputPrecheck.load_checkpoint()
        assert isinstance(params, dict)
        assert isinstance(wrapper, MFA2TreeWrapper)

    def test_save_checkpoint(self, tmp_path):
        OutputPrecheck.setup(folder=tmp_path)
        params = {
            "method": "ft",
            "top_n_toverr": 18,
            "concat": False,
            "partition": None,
            "inputs": "fakechecksum",
        }
        wrapper = MFA2TreeWrapper(["tests/data/pep_align/10at10240.aa.mfa"], seqtype="pep")
        OutputPrecheck.save_checkpoint(params, wrapper)
        a, b = OutputPrecheck.load_checkpoint()
        assert len(a) == 6
        assert b == wrapper

    @pytest.mark.parametrize(
        "method, concat, partition, figure",
        (
            ("ft", True, "seq", False),
            ("iqtree", False, None, False),
            ("raxml", True, "seq+codon", True),
        ),
    )
    def test_output_results(self, copy_phylotree_ckp: Path, method: str, concat: bool, partition: str, figure: bool):
        tmp_path = copy_phylotree_ckp
        OutputPrecheck.setup(folder=tmp_path, method=method, concat=concat, partition=partition, figure=figure)
        _, wrapper = OutputPrecheck.load_checkpoint()
        OutputPrecheck.output_results(wrapper[0].tree)

        with open(tmp_path / "final_tree.nw") as f:
            content = f.read()
        strategy = "concatenate" if concat else "consensus"
        assert f"Final tree is built using {_config.avail_tree_methods[method]} with {strategy} strategy" in content

        if concat and partition:
            assert f"Partition is enabled using {partition} mode" in content

        if figure:
            assert (tmp_path / "final_tree.png").is_file()


@pytest.mark.parametrize("input_dir, exp_seqtype", ((Path("tests/data/pep_align"), "pep"), (Path("tests/data/cds_align"), "cds")))
def test_determine_samples_and_seqtype(input_dir, exp_seqtype):
    samples, seqtype = determine_samples_and_seqtype(input_dir)

    assert len(samples) == 5
    assert seqtype == exp_seqtype


def test_determine_samples_and_seqtype_checkpoint_not_found(tmp_path, caplog: pytest.LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        samples, seqtype = determine_samples_and_seqtype(tmp_path)

    assert samples is None
    assert seqtype is None
    assert "Determine samples and seqtype from the given inputs" in caplog.text
