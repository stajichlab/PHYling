"""IQ-Tree utilities"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Literal, overload

from ..libphyling import SeqTypes, TreeMethods
from ..libphyling._utils import check_binary
from ._base import BinaryWrapper, TreeToolWrapper

IQTREE_BIN = check_binary(
    TreeMethods.IQTREE.method, TreeMethods.IQTREE.bins, "bioconda::iqtree", "https://github.com/iqtree/iqtree2"
)


def partition_nexus_to_raxml(file: str | Path, output: str | Path):
    with open(file) as f:
        sections = f.read().strip("\n").split(";\n")
        info = {}
        for charset in sections[1:-2]:
            charset, pos = re.sub("charset", "", charset).strip().split(" = ")
            info[charset] = [pos]
        for model in re.sub(r".*\n", "", sections[-2], count=1).split(",\n"):
            m, charset = model.strip().split(": ")
            m = re.sub(r"\+G4", "+G", m)
            info[charset].append(m)
    with open(output, "w") as f:
        for k, v in info.items():
            f.write(f"{v[1]}, {k} = {v[0]}\n")


class ModelFinder(BinaryWrapper):
    @overload
    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        partition_file: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        method: Literal["ft", "raxml", "iqtree"] = "iqtree",
        threads: int = 1,
        add_args: tuple | list | None = None,
    ): ...

    @overload
    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        method: Literal["ft", "raxml", "iqtree"] = "iqtree",
        threads: int = 1,
        add_args: tuple | list | None = None,
    ): ...

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        partition_file: str | Path | None = None,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        method: Literal["ft", "raxml", "iqtree"] = "iqtree",
        threads: int = 1,
        add_args: tuple | list | None = None,
    ):
        super().__init__(
            "ModelFinder", file, output, partition_file, seqtype=seqtype, method=method, add_args=add_args, threads=threads
        )

    @overload
    def __call__(self, *, verbose=False) -> Path: ...
    @overload
    def __call__(self, *, verbose=False) -> str: ...

    def __call__(self, *, verbose=False) -> str | Path:
        super().__call__(verbose=verbose)
        if self._target.suffix == ".nex":  # partitioning analysis
            partition_nexus_to_raxml(self._target, self._target.with_suffix(""))
            return self._target.with_suffix("")
        else:
            with gzip.open(self._target, "rt") as f:
                for line in f.read().strip("\n").split("\n"):
                    best_model_prefix = "best_model_BIC: "
                    if line.startswith(best_model_prefix):
                        best_model = line.lstrip(best_model_prefix)
                        if self._method == "ft":
                            best_model = re.sub(r"\+I(?=\+)", "", re.sub(r"\+F(?=\+)", "", best_model))
            return best_model

    def _params_check(
        self,
        partition_file: str | Path,
        seqtype: Literal["dna", "pep", "AUTO"],
        method: Literal["ft", "raxml", "iqtree"],
        **kwargs,
    ):
        if partition_file:
            if not isinstance(partition_file, (str, Path)):
                raise TypeError(f"Argument partition only accepts str or Path. Got {type(partition_file)}.")
            if method == "ft":
                raise ValueError(f"Partitioning analysis is not allowed when using {TreeMethods.FT.method}.")
            partition_file = Path(partition_file)
        if seqtype == SeqTypes.DNA:
            seqtype = "DNA"
            mset = ",".join(TreeMethods[method.upper()].dna_model)
        elif seqtype == SeqTypes.PEP:
            seqtype = "AA"
            mset = ",".join(TreeMethods[method.upper()].pep_model)
        else:
            seqtype = None
            mset = None
        self._method = method
        return super()._params_check(partition_file, seqtype=seqtype, mset=mset, **kwargs)

    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        partition_file: Path | None,
        *,
        seqtype: Literal["dna", "pep"] | None,
        mset: str | None,
        threads: int,
    ):
        self._cmd = [
            IQTREE_BIN,
            "-s",
            str(file.absolute()),
            "--prefix",
            str(output.absolute()),
            "-T",
            "AUTO",
            "--threads-max",
            str(threads),
            "-m",
            "TESTONLY",
        ]
        if seqtype:
            self._cmd.extend(["--seqtype", seqtype])
        if mset:
            self._cmd.extend(["--mset", mset])
        if partition_file:
            self._cmd.extend(["-p", str(partition_file)])
            self._target = Path(f"{output}.best_scheme.nex")
        else:
            self._target = Path(f"{output}.model.gz")


class Iqtree(TreeToolWrapper):
    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        model: str = "AUTO",
        threads: int = 1,
        add_args: tuple | list | None = None,
    ):
        super().__init__(
            TreeMethods.IQTREE.method, file, output, seqtype=seqtype, model=model, add_args=add_args, threads=threads
        )

    def _construct_cmd(self, file: Path, output: Path, *, seqtype: Literal["DNA", "AA"] | None, model: str, threads: int):
        self._cmd = [
            IQTREE_BIN,
            "-s",
            str(file.absolute()),
            "--prefix",
            str(output.absolute()),
            "-T",
            "AUTO",
            "--threads-max",
            str(threads),
        ]
        if seqtype:
            self._cmd.extend(["--seqtype", seqtype])
        if Path(model).is_file():
            self._cmd.extend(["-p", str(model)])
        else:
            self._cmd.extend(["-m", str(model)])

        self._target = Path(f"{output}.treefile")

    def _update_model(self):
        model_file = self._target.with_suffix(".best_model.nex")

        if model_file.is_file():
            self._model = model_file
        else:
            with open(self._target.with_suffix(".iqtree")) as f:
                self._model = re.search(r"alisim simulated_MSA .* (\-m) \"(.*)\" ", f.read())[2]


class UFBoot(TreeToolWrapper):
    def __init__(
        self,
        file: str | Path,
        tree: str | Path,
        output: str | Path,
        *,
        model: str = "AUTO",
        bs: int = 1000,
        threads: int = 1,
        add_args: tuple | list | None = None,
    ):
        super().__init__(
            "Ultrafast Bootstrap", file, output, tree, seqtype="AUTO", model=model, add_args=add_args, threads=threads, bs=bs
        )

    def _params_check(self, tree: str | Path, *, seqtype, **kwargs):
        tree = Path(tree)
        if not tree.exists():
            raise FileNotFoundError(f"{tree}")
        if not tree.is_file():
            raise RuntimeError(f"{tree} is not a file.")
        return super()._params_check(tree, seqtype=seqtype, **kwargs)

    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        tree: Path,
        *,
        seqtype: None,
        model: str,
        bs: int,
        threads: int,
    ):
        self._cmd = [
            IQTREE_BIN,
            "-s",
            str(file.absolute()),
            "--prefix",
            str(output.absolute()),
            "-t",
            str(tree.absolute()),
            "-B",
            str(bs),
            "-T",
            "AUTO",
            "--threads-max",
            str(threads),
        ]
        if Path(model).is_file():
            self._cmd.extend(["-p", str(model)])
        else:
            self._cmd.extend(["-m", str(model)])

        self._target = Path(f"{output}.treefile")


class Concordance(TreeToolWrapper):
    def __init__(
        self,
        file: str | Path,
        tree: str | Path,
        output: str | Path,
        *,
        model: str = "AUTO",
        scfl: int = 100,
        threads: int = 1,
        add_args: tuple | list | None = None,
    ):
        super().__init__(
            "Concordance factor calculation",
            file,
            output,
            tree,
            seqtype="AUTO",
            model=model,
            add_args=add_args,
            threads=threads,
            scfl=scfl,
        )

    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        tree: Path,
        *,
        seqtype: None,
        model: str,
        scfl: int,
        threads: int,
    ):
        self._cmd = [
            IQTREE_BIN,
            "-s",
            str(file.absolute()),
            "--prefix",
            str(output.absolute()),
            "-te",
            str(tree.absolute()),
            "--scfl",
            "100",
            "-T",
            "AUTO",
            "--threads-max",
            str(threads),
        ]
        if Path(model).is_file():
            self._cmd.extend(["-p", str(model)])
        else:
            self._cmd.extend(["-m", str(model)])

        self._target = Path(f"{output}.cf.tree")
