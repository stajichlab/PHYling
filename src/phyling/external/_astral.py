"""Astral utilities"""

from __future__ import annotations

from pathlib import Path

from ..libphyling._utils import check_binary
from ._abc import BinaryWrapper

ASTRAL_BIN = check_binary("ASTRAL", ("astral",), "bioconda::aster", "https://github.com/chaoszhang/ASTER")


class Astral(BinaryWrapper):
    """Compute a consensus tree using ASTRAL.

    Returns:
        Tree: The consensus tree.

    Raises:
        BinaryNotFoundError: If ASTRAL is not installed.
        RuntimeError: If ASTRAL fails.
    """

    _prog: str = "ASTRAL"
    _cmd_log = "stderr"

    def __init__(self, file: str | Path, output: str | Path, *, add_args: list | tuple | None = None, threads: int = 1):
        super().__init__(file, output, add_args=add_args, threads=threads)

    def _construct_cmd(self, file: Path, output: Path, *, threads: int):
        self._cmd = [ASTRAL_BIN, "--output", str(output), "--thread", str(threads), str(file)]
