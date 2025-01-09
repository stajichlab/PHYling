from __future__ import annotations

import logging
from argparse import ArgumentParser
from itertools import product
from os import cpu_count
from pathlib import Path

import pytest

from phyling.libphyling import ALIGN_METHODS, TreeMethods, TreeOutputFiles
from phyling.pipeline import align, download, filter, tree


@pytest.fixture
def parser() -> ArgumentParser:
    """Fixture to create a new parser for each test."""
    return ArgumentParser(add_help=False)


def is_subset(d1: dict, d2: dict):
    return all(k in d1 and d1[k] == v for k, v in d2.items())


@pytest.mark.usefixtures("hide_metadata_wo_download")
class TestDownload:
    def test_menu(self, parser: ArgumentParser):
        download.menu(parser)
        true_param_pairs = vars(parser.parse_args(["poxviridae_odb10"]))
        expected_param_pairs = {"markerset": "poxviridae_odb10"}
        assert is_subset(true_param_pairs, expected_param_pairs)

    def test_download_list(self):
        download.download("list")

    def test_download_with_valid_name(self, capsys: pytest.CaptureFixture):
        with capsys.disabled():
            download.download("poxviridae_odb10")
        download.download("list")
        captured: str = capsys.readouterr().out
        captured = captured.split("Datasets available on local:")
        assert "poxviridae_odb10" in captured[1]

    def test_download_with_already_up_to_date(self, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.INFO):
            download.download("poxviridae_odb10")
        assert "Markerset already exists and update to date." in caplog.text

    def test_download_with_invalid_name(self):
        with pytest.raises(RuntimeError) as excinfo:
            download.download("InvalidName")
        assert "Markerset not available" in str(excinfo.value)


class TestAlign:
    inputs = (
        "tests/data/pep",
        "tests/data/cds",
        [
            "tests/data/pep/Anomala_cuprea_entomopoxvirus.faa.gz",
            "tests/data/pep/Canarypox_virus.faa.gz",
            "tests/data/pep/Cowpox_virus.faa.gz",
            "tests/data/pep/Goatpox_virus.faa.gz",
            "tests/data/pep/Variola_virus.faa.gz",
        ],
        [
            "tests/data/cds/Anomala_cuprea_entomopoxvirus.fna.gz",
            "tests/data/cds/Canarypox_virus.fna.gz",
            "tests/data/cds/Cowpox_virus.fna.gz",
            "tests/data/cds/Goatpox_virus.fna.gz",
            "tests/data/cds/Variola_virus.fna.gz",
        ],
    )
    nontrim = (True, False)

    def test_menu(self, parser: ArgumentParser):
        align.menu(parser)
        true_param_pairs = vars(
            parser.parse_args(
                [
                    "--input_dir",
                    "I",
                    "--output",
                    "O",
                    "--markerset",
                    "M",
                    "--evalue",
                    "1e-10",
                    "--method",
                    "hmmalign",
                    "--non_trim",
                    "--threads",
                    "3",
                ]
            )
        )
        expected_param_pairs = {
            "inputs": Path("I"),
            "output": Path("O"),
            "markerset": Path("M"),
            "evalue": 1e-10,
            "method": "hmmalign",
            "non_trim": True,
            "threads": 3,
        }
        assert is_subset(true_param_pairs, expected_param_pairs)

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method, nontrim", list(product(inputs, ALIGN_METHODS, nontrim)))
    def test_align_inputs(self, inputs: str | list[str], tmp_path: Path, method: str, nontrim: bool):
        align.align(
            inputs,
            tmp_path,
            markerset="tests/database/poxviridae_odb10/hmms",
            method=method,
            non_trim=nontrim,
            threads=1,
        )
        filenames = [file.name for file in tmp_path.iterdir()]
        first_file = inputs if isinstance(inputs, str) else inputs[0]

        if first_file.startswith("tests/data/pep"):
            assert filenames == [file.name for file in Path("tests/data/pep_align").iterdir()]
        else:
            assert filenames == [file.name for file in Path("tests/data/cds_align").iterdir()]

    def test_align_too_few_inputs(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            align.align(self.inputs[2][2:], tmp_path, markerset="tests/database/poxviridae_odb10/hmms")
        assert "Requires at least 4 samples" in str(excinfo.value)

    def test_align_invalid_markerset(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError) as excinfo:
            align.align(self.inputs[0], tmp_path, markerset="InvalidName")
        assert "Markerset folder does not exist" in str(excinfo.value)

    def test_align_invalid_evalue(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            align.align(self.inputs[0], tmp_path, markerset="tests/database/poxviridae_odb10/hmms", evalue=1)
        assert "Invalid evalue" in str(excinfo.value)

    def test_align_invalid_method(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            align.align(self.inputs[0], tmp_path, markerset="tests/database/poxviridae_odb10/hmms", method="InvalidName")
        assert "Invalid method" in str(excinfo.value)


class TestFilter:
    inputs_pep = (
        "tests/data/pep_align",
        [
            "tests/data/pep_align/10at10240.aa.mfa",
            "tests/data/pep_align/14at10240.aa.mfa",
            "tests/data/pep_align/22at10240.aa.mfa",
            "tests/data/pep_align/43at10240.aa.mfa",
            "tests/data/pep_align/48at10240.aa.mfa",
            "tests/data/pep_align/53at10240.aa.mfa",
        ],
    )
    inputs_cds = (
        "tests/data/cds_align",
        [
            "tests/data/cds_align/10at10240.cds.mfa",
            "tests/data/cds_align/14at10240.cds.mfa",
            "tests/data/cds_align/22at10240.cds.mfa",
            "tests/data/cds_align/43at10240.cds.mfa",
            "tests/data/cds_align/48at10240.cds.mfa",
            "tests/data/cds_align/53at10240.cds.mfa",
        ],
    )
    invalid_top_n = (1, 7)

    def test_menu(self, parser: ArgumentParser):
        filter.menu(parser)
        true_param_pairs = vars(
            parser.parse_args(["--input_dir", "I", "--output", "O", "--top_n_toverr", "10", "--threads", "3"])
        )
        expected_param_pairs = {
            "inputs": Path("I"),
            "output": Path("O"),
            "top_n_toverr": 10,
            "threads": 3,
        }
        assert is_subset(true_param_pairs, expected_param_pairs)

    @pytest.mark.parametrize("inputs", inputs_pep)
    def test_filter_pep(self, inputs: str | list[str], tmp_path: Path):
        filter.filter(inputs, tmp_path, top_n_toverr=5, threads=cpu_count())
        treeness = tmp_path / TreeOutputFiles.TREENESS
        msas_dir = tmp_path / TreeOutputFiles.MSAS_DIR
        assert treeness.exists()
        assert treeness.is_file()
        assert msas_dir.exists()
        assert msas_dir.is_dir()

    @pytest.mark.parametrize("inputs", inputs_cds)
    def test_filter_cds(self, inputs: str | list[str], tmp_path: Path):
        filter.filter(inputs, tmp_path, top_n_toverr=5, threads=cpu_count())
        treeness = tmp_path / TreeOutputFiles.TREENESS
        msas_dir = tmp_path / TreeOutputFiles.MSAS_DIR
        assert treeness.exists()
        assert treeness.is_file()
        assert msas_dir.exists()
        assert msas_dir.is_dir()

    @pytest.mark.parametrize("input", list((inputs_pep[1][:2], inputs_cds[1][:2])))
    def test_filter_too_few_input(self, input: str, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            filter.filter(input, tmp_path, top_n_toverr=5)
        assert "Fewer than 3 inputs" in str(excinfo.value)

    @pytest.mark.parametrize("input, top_n", list(product((inputs_pep[1], inputs_pep[1][:3]), invalid_top_n)))
    def test_filter_invalid_top_n(self, input: str, tmp_path: Path, top_n: int):
        with pytest.raises(ValueError) as excinfo:
            filter.filter(input, tmp_path, top_n_toverr=top_n)
        assert "Argument top_n_toverr out of range" in str(excinfo.value)

    @pytest.mark.parametrize("input", (inputs_pep[1],))
    def test_filter_top_n_equal_to_inputs(self, input: str, tmp_path: Path):
        with pytest.raises(SystemExit) as excinfo:
            filter.filter(input, tmp_path, top_n_toverr=6)
        assert "Argument top_n_toverr is equal to the number of inputs" in str(excinfo.value)


class TestTree:
    inputs_pep = (
        "tests/data/pep_align",
        [
            "tests/data/pep_align/10at10240.aa.mfa",
            "tests/data/pep_align/14at10240.aa.mfa",
            "tests/data/pep_align/43at10240.aa.mfa",
        ],
    )
    inputs_cds = (
        "tests/data/cds_align",
        [
            "tests/data/cds_align/10at10240.cds.mfa",
            "tests/data/cds_align/14at10240.cds.mfa",
            "tests/data/cds_align/22at10240.cds.mfa",
        ],
    )
    method = tuple(m.name.lower() for m in TreeMethods)
    concat = (False, True)
    partition = (False, True)
    figure = (False, True)

    def test_menu(self, parser: ArgumentParser):
        tree.menu(parser)
        true_param_pairs = vars(
            parser.parse_args(
                [
                    "--input_dir",
                    "I",
                    "--output",
                    "O",
                    "--method",
                    "iqtree",
                    "--bootstrap",
                    "100",
                    "--concat",
                    "--partition",
                    "--threads",
                    "3",
                ]
            )
        )
        expected_param_pairs = {
            "inputs": Path("I"),
            "output": Path("O"),
            "method": "iqtree",
            "bootstrap": 100,
            "concat": True,
            "partition": True,
            "threads": 3,
        }
        assert is_subset(true_param_pairs, expected_param_pairs)

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method, concat", list(product(inputs_pep, method, concat)))
    def test_tree_pep(self, inputs: str | list[str], tmp_path: Path, method: str, concat: bool):
        tree.tree(inputs, tmp_path, method=method, bs=10, concat=concat, threads=1)
        treefile = tmp_path / TreeOutputFiles.TREE_NW
        assert treefile.exists()
        assert treefile.is_file()

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method, concat", list(product(inputs_cds, method, concat)))
    def test_tree_cds(self, inputs: str | list[str], tmp_path: Path, method: str, concat: bool):
        tree.tree(inputs, tmp_path, method=method, bs=10, concat=concat, threads=1)
        treefile = tmp_path / TreeOutputFiles.TREE_NW
        assert treefile.exists()
        assert treefile.is_file()

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method", list(product((inputs_pep[1],), method[1:])))
    def test_tree_pep_partition(self, inputs: str | list[str], tmp_path: Path, method: str):
        tree.tree(inputs, tmp_path, method=method, bs=10, concat=True, partition=True, threads=1)
        treefile = tmp_path / TreeOutputFiles.TREE_NW
        assert treefile.exists()
        assert treefile.is_file()

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method", list(product((inputs_cds[1],), method[1:])))
    def test_tree_cds_partition(self, inputs: str | list[str], tmp_path: Path, method: str):
        tree.tree(inputs, tmp_path, method=method, bs=10, concat=True, partition=True, threads=1)
        treefile = tmp_path / TreeOutputFiles.TREE_NW
        assert treefile.exists()
        assert treefile.is_file()

    @pytest.mark.parametrize("inputs, method", list(product((inputs_pep[1],), method[1:])))
    def test_partition_no_concat(self, inputs: str | list[str], tmp_path: Path, method: str):
        with pytest.raises(ValueError) as excinfo:
            tree.tree(inputs, tmp_path, method=method, bs=10, concat=False, partition=True)
        assert "Partition is not allowed in consensus mode" in str(excinfo.value)

    def test_partition_invalid_method(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            tree.tree(self.inputs_pep[1], tmp_path, bs=10, concat=True, partition=True)
        assert f"Partition is not allowed with {TreeMethods.FT.method}" in str(excinfo.value)

    @pytest.mark.parametrize("input, method", list(product((inputs_pep[1][0], inputs_cds[1][0]), method)))
    def test_tree_single_input(self, input: list, tmp_path: Path, method: str):
        with pytest.raises(ValueError) as excinfo:
            tree.tree(input, tmp_path, bs=10, method=method)
        assert "Found only 1 MSA fasta" in str(excinfo.value)

    @pytest.mark.parametrize("inputs", inputs_pep)
    def test_tree_figure(self, inputs, tmp_path: Path):
        tree.tree(inputs, tmp_path, bs=10, figure=True)
        assert tmp_path / TreeOutputFiles.TREE_IMG
