import subprocess
from pathlib import Path

from .. import logger
from ..libphyling._utils import check_binary
from ._base import BinaryWrapper

ASTRAL_BIN = check_binary("ASTRAL", ("astral",), "bioconda::aster", "https://github.com/chaoszhang/ASTER")


class Astral(BinaryWrapper):
    """Compute a consensus tree using ASTRAL.

    Returns:
        Tree: The consensus tree.

    Raises:
        BinaryNotFoundError: If ASTRAL is not installed.
        RuntimeError: If ASTRAL fails.
    """

    def __init__(self, file: str | Path, output: str | Path, *, add_args: list | tuple | None = None, threads: int = 1):
        super().__init__("ASTRAL", file, output, add_args=add_args, threads=threads)

    def __call__(self, *, verbose: bool = False) -> Path:
        try:
            result = subprocess.run(self._cmd, capture_output=True, check=True, text=True)
            if verbose:
                logger.debug("%s", result.stderr)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"{self._prog} failed with cmd: {self.cmd}:\n{e.stderr}")
        return self._target

    def _construct_cmd(self, file: Path, output: Path, *, threads: int):
        self._cmd = [ASTRAL_BIN, "--output", str(output), "--thread", str(threads), str(file)]
        self._target = output
