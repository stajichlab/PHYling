"""FastTree utilities"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Literal

from .. import logger
from ..libphyling import TreeMethods
from ..libphyling._utils import check_binary
from ._base import TreeToolWrapper

FASTTREE_BIN = check_binary(
    TreeMethods.FT.method, TreeMethods.FT.bins, "bioconda::veryfasttree", "https://github.com/citiususc/veryfasttree"
)


class FastTree(TreeToolWrapper):
    """Runs FastTree to build a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree (MFA2Tree): The MFA2Tree object containing multiple sequence alignment data.
        capture_cmd (bool): If True, returns the FastTree command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple containing the Tree object and the command string.
    """

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep"],
        model: str = "AUTO",
        threads: int = 1,
        add_args: tuple | list | None = None,
    ):
        super().__init__(TreeMethods.FT.method, file, output, seqtype=seqtype, model=model, add_args=add_args, threads=threads)

    def __call__(self, *, verbose: bool = False) -> Path:
        try:
            result = subprocess.run(self._cmd, capture_output=True, check=True, text=True)
            if verbose:
                logger.debug("%s", result.stderr)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"{self._prog} failed with cmd: {self.cmd}:\n{e.stderr}")
        return self._target

    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        *,
        seqtype: Literal["DNA", "AA"],
        model: str = "AUTO",
        threads: int = 1,
    ):
        self._cmd = [
            FASTTREE_BIN,
            "-nosupport",
            "-out",
            str(output),
            "-threads",
            str(threads),
            str(file),
        ]
        model = model.split("+")
        if seqtype == "DNA":
            self._cmd.insert(1, "-nt")
            if model[0].upper() not in TreeMethods.FT.dna_model:
                raise ValueError(f"Model {model[0]} is not supported in {TreeMethods.FT.method} with {seqtype} alignments.")
            if model[0].lower() in TreeMethods.FT.dna_model[1:]:
                self._cmd.insert(3, f"--{model[0].lower()}")
        else:
            if model[0].upper() not in TreeMethods.FT.pep_model:
                raise ValueError(f"Model {model[0]} is not supported in {TreeMethods.FT.method} with {seqtype} alignments.")
            if model[0].lower() in TreeMethods.FT.pep_model[1:]:
                self._cmd.insert(2, f"--{model[0].lower()}")

        for model_param in model[1:]:
            if model_param == "G":
                self._cmd.insert(-1, "-gamma")

        self._target = output
