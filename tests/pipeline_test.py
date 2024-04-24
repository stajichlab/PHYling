import logging
from itertools import product
from pathlib import Path

import pytest

from phyling.pipeline import align, download, tree


@pytest.mark.usefixtures("hide_metadata_wo_download")
class TestDownload:
    def test_download_list(self):
        download("list")

    def test_download_with_valid_name(self, capsys: pytest.CaptureFixture):
        with capsys.disabled():
            download("poxviridae_odb10")
        download("list")
        captured: str = capsys.readouterr().out
        captured = captured.split("Datasets available on local:")
        assert "poxviridae_odb10" in captured[1]

    def test_download_with_already_up_to_date(self, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.INFO):
            download("poxviridae_odb10")
        assert "Markerset already exists and update to date." in caplog.text

    def test_download_with_invalid_name(self, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.ERROR):
            with pytest.raises(SystemExit):
                download("InvalidName")

        assert any(record.levelname == "ERROR" for record in caplog.records)
        assert "Markerset InvalidName not available" in caplog.text


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
    method = ("hmmalign", "muscle")
    nontrim = (True, False)

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method, nontrim", list(product(inputs, method, nontrim)))
    def test_align_inputs(self, inputs: str, tmp_path: Path, method: str, nontrim: bool):
        align(inputs, tmp_path, markerset="tests/database/poxviridae_odb10/hmms", method=method, non_trim=nontrim)

    def test_align_markerset_not_avail(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        with pytest.raises(SystemExit):
            align(self.inputs[0], tmp_path, markerset="InvalidName")

        assert any(record.levelname == "ERROR" for record in caplog.records)
        assert "Markerset folder does not exist" in caplog.text

    def test_align_method_not_avail(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        with pytest.raises(SystemExit):
            align(self.inputs[0], tmp_path, markerset="tests/database/poxviridae_odb10/hmms", method="InvalidName")

        assert any(record.levelname == "ERROR" for record in caplog.records)
        assert "Invalid method" in caplog.text


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
    method = ("ft", "raxml", "iqtree")
    top_n = (0, 4, 50)
    concat = (False, True)
    partition = (None, "seq", "codon", "seq+codon")
    figure = (False, True)

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, top_n, concat, partition", list(product(inputs_pep, top_n, concat, partition[:2])))
    def test_tree_pep(
        self, inputs: str | list[str], tmp_path: Path, top_n: int, concat: bool, partition: str, caplog: pytest.LogCaptureFixture
    ):
        with caplog.at_level(logging.WARNING):
            tree(inputs, tmp_path, top_n_toverr=top_n, concat=concat, partition=partition)

        msg = "Partition is forced to be disabled since it only works when using raxml and iqtree."
        if partition:
            assert msg in caplog.text

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, top_n, concat, partition", list(product(inputs_cds, top_n, concat, partition)))
    def test_tree_cds(
        self, inputs: str | list[str], tmp_path: Path, top_n: int, concat: bool, partition: str, caplog: pytest.LogCaptureFixture
    ):
        with caplog.at_level(logging.WARNING):
            tree(inputs, tmp_path, top_n_toverr=top_n, concat=concat, partition=partition)

        msg = "Partition is forced to be disabled since it only works when using raxml and iqtree."
        if partition:
            assert msg in caplog.text

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "inputs, method, concat, partition",
        list(product((inputs_pep[1], inputs_cds[1]), method[1:], concat, partition[:2])),
    )
    def test_tree_raxml_iqtree(
        self, inputs: list, tmp_path: Path, method: str, concat: bool, partition: str, caplog: pytest.LogCaptureFixture
    ):
        with caplog.at_level(logging.WARNING):
            tree(inputs, tmp_path, method=method, concat=concat, partition=partition)

        msg = 'Peptide sequence detected. Force to use "seq" mode for partitioning.'
        if inputs in self.inputs_pep and concat and partition in self.partition[2:]:
            assert msg in caplog.text
        else:
            assert msg not in caplog.text

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "input, method",
        list(product((inputs_pep[1][0], inputs_cds[1][0]), method)),
    )
    def test_tree_single_input(self, input: list, tmp_path: Path, method: str, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.INFO):
            tree(input, tmp_path, method=method)

        assert "Only one MSA is selected. Report the tree directly." in caplog.text

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "inputs, method, concat",
        list(product((inputs_pep[1], inputs_cds[1]), method[1:], concat)),
    )
    def test_tree_top_1(self, inputs: list, tmp_path: Path, method: str, concat: bool, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.ERROR):
            with pytest.raises(SystemExit):
                tree(inputs, tmp_path, method=method, top_n_toverr=1, concat=concat)

        assert "cannot be run since only single mfa found after filtering by toverr." in caplog.text

    @pytest.mark.parametrize("inputs", inputs_pep)
    def test_tree_figure(self, inputs, tmp_path: Path):
        tree(inputs, tmp_path, figure=True)
        assert tmp_path / "inal_tree.png"
