import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Literal

from .. import logger
from ..libphyling import SeqTypes


class BinaryWrapper(ABC):
    def __init__(
        self, prog: str, file: str | Path, output: str | Path, *args, add_args: list | tuple | None = None, **kwargs
    ) -> None:
        self._prog = prog
        file = Path(file)
        if not file.exists():
            raise FileNotFoundError(f"{file}")
        output = Path(output)
        add_args = add_args if add_args else []
        if not isinstance(add_args, (list, tuple)):
            raise TypeError(f"Argument add_args only accepts list or tuple. Got {type(self.add_args)}.")
        args, kwargs = self._params_check(*args, **kwargs)
        self._target: Path = None
        self._construct_cmd(file, output, *args, **kwargs)
        self._cmd: list[str]
        self._cmd.extend(add_args)

    def __call__(self, *, verbose: bool = False) -> Path:
        """Execute the command.

        Returns:
            Path: The path of the self._target.
        """
        if self._target:
            self._target.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(self.cmd)
        try:
            result = subprocess.run(self._cmd, capture_output=True, check=True, text=True)
            if verbose:
                logger.debug("%s", result.stdout)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"{self._prog} failed with cmd: {self.cmd}\n{e.stderr}")
        if not self._target:
            return result.stdout
        else:
            return self._target

    @property
    def cmd(self) -> str:
        return " ".join(self._cmd)

    def _params_check(self, *args, **kwargs) -> tuple[tuple, dict]:
        return args, kwargs

    @abstractmethod
    def _construct_cmd(self, file: Path, output: Path, *args, **kwargs) -> None: ...


class TreeToolWrapper(BinaryWrapper):
    def __init__(
        self,
        prog: str,
        file: str | Path,
        output: str | Path,
        *args,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        model: str = "AUTO",
        add_args: list | tuple | None = None,
        threads: int = 1,
        **kwargs,
    ) -> None:
        super().__init__(prog, file, output, *args, seqtype=seqtype, model=model, add_args=add_args, threads=threads, **kwargs)
        self._model: str = model

    def __call__(self, *, verbose=False):
        """Execute the phylogeny inference command.

        Returns:
            Path: The path of the newick tree file.
        """
        result = super().__call__(verbose=verbose)
        self._update_model()
        return result

    @property
    def model(self):
        return self._model

    def _params_check(self, *args, seqtype: str, **kwargs) -> tuple[tuple, dict]:
        if seqtype == SeqTypes.DNA:
            seqtype = "DNA"
        elif seqtype == SeqTypes.PEP:
            seqtype = "AA"
        else:
            seqtype = None
        return super()._params_check(*args, seqtype=seqtype, **kwargs)

    @abstractmethod
    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        *args,
        seqtype: Literal["DNA", "AA"] | None,
        model: str,
        threads: int,
        **kwargs,
    ) -> None: ...

    def _update_model(self):
        pass
