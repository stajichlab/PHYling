from pathlib import Path

from ..libphyling._utils import check_binary
from ._base import BinaryWrapper

MUSCLE_BIN = check_binary("Muscle", ("muscle",), "bioconda::muscle", "https://github.com/rcedgar/muscle")


class Muscle(BinaryWrapper):
    def __init__(self, file: str | Path, output: str | Path, threads: int = 1, add_args: tuple | list | None = None):
        super().__init__("Muscle", file, output, add_args=add_args, threads=threads)

    def _construct_cmd(self, file: Path, output: Path, *, threads: int = 1):
        self._cmd = [
            MUSCLE_BIN,
            "-align",
            str(file),
            "-output",
            str(output),
            "-threads",
            str(threads),
        ]

        self._target = output
