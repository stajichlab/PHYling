import logging
import pickle
import shutil
import textwrap
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, Iterable, Iterator, TypeVar

import phyling.config as config

T = TypeVar("T")
SFW = TypeVar("SFW", bound="SeqFileWrapperABC")


class SeqFileWrapperABC(ABC):
    """An abstract class for sequence files."""

    def __init__(self, file: str | Path, seqtype: str | None = None) -> None:
        """Initialize the object and check sequence type."""
        self._file = Path(file)
        if seqtype:
            self._set_seqtype(seqtype)

    def __repr__(self) -> str:
        """Return the string representation."""
        return f"{self.__class__.__qualname__}({self.name}" + ", %s)"

    @abstractmethod
    def __eq__(self, other: SFW) -> bool:
        """Return true if the two objects have the same value."""
        if not isinstance(other, self.__class__):
            raise TypeError(f"Can only compare to {self.__class__} object but got {type(other)}")

    @abstractmethod
    def __lt__(self, other: SFW) -> bool:
        """Return true if the current object is smaller than/prior to another object."""
        ...

    @property
    def path(self) -> Path:
        """Return the pathlib object of the file."""
        return self._file

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of the entry."""
        ...

    @property
    def seqtype(self) -> str:
        """Return whether the sequences are peptide or DNA."""
        return self._seqtype

    def _set_seqtype(self, seqtype: str):
        """Set the self._seqtype."""
        supported_seqtype = (config.seqtype_pep, config.seqtype_cds)
        if seqtype not in supported_seqtype:
            raise KeyError(f'Argument "seqtype" not falling in {supported_seqtype}')
        self._seqtype = seqtype


class DataListABC(ABC, Generic[T]):
    """A list-like abstract class which have partial list methods."""

    @abstractmethod
    def __init__(self, files: Iterable[str | Path | T]) -> None:
        """Initialize the object and store data into a list."""
        self.data: list[T] = []
        if not isinstance(files, Iterable):
            raise TypeError(f"{self.__class__.__qualname__} is not iterable")

    def __repr__(self) -> str:
        """Return the object representation when being called."""
        if len(repr(self.data)) < 40:
            data = self.data
        else:
            data = [str(data) for data in self.data]
            if len(self) > 20:
                data = data[:5] + ["..."] + data[-5:]
            data = ",\n".join([textwrap.indent(line, prefix=" " * 4) for line in data])
            data = f"[\n{data}]"
        return self.__class__.__qualname__ + "%s" + data

    def __iter__(self) -> Iterator[T]:
        """Make the object iterable."""
        return iter(self.data)

    def __eq__(self, other) -> bool:
        """Return true if the checksums of two objects are the same."""
        if not isinstance(other, self.__class__):
            raise TypeError(f"Can only compare to {self.__class__.__qualname__} object but got {type(other)}.")
        return self.checksum == other.checksum

    def __contains__(self, item: T) -> None:
        """Return true if item is contained in the object."""
        return item in self.data

    def __len__(self) -> int:
        """Return the number of data contain in this object."""
        return len(self.data)

    def __getitem__(self, idx: int | slice) -> T | list[T]:
        """Get the data by index or slice."""
        if isinstance(idx, slice):
            return self.__class__(self.data[idx])
        else:
            return self.data[idx]

    @property
    @abstractmethod
    def checksum(self) -> str:
        """Compute the checksum of the stored data."""
        ...

    def sort(self, /, *args, **kwds) -> None:
        """Sort the data collection."""
        self.data.sort(*args, **kwds)

    def _names(self) -> tuple[str]:
        """Get the tuple of the names of the stored objects."""
        return tuple(data.name for data in self.data)


class OutputPrecheck(ABC, Generic[T]):
    """A abstract class that provides features for input/output precheck and checkpoint loading/saving."""

    @abstractmethod
    def __init__(self, output: str | Path, param_keys: tuple[str], data: T = None) -> None:
        """Initialize the object with Nonetype slots for parameters and data if not receiving params and samplelist."""
        self._output = Path(output)
        self._params = {param: None for param in param_keys}
        self._data = data
        self._ckp: Path  # Need to be specified in the subclass.
        ...

    def __str__(self) -> str:
        """Return the string representation of the object."""
        return str(self._output)

    @property
    def path(self) -> Path:
        """Return the pathlib object of the file."""
        return self._output

    @property
    def checkpoint(self) -> Path:
        """Return the pathlib object of the checkpoint."""
        return self._ckp

    @property
    def params(self) -> dict:
        """Return the parameters of the current run."""
        return self._params

    @params.setter
    def params(self, cur_params: dict) -> None:
        """Update the current parameters by the given dictionary."""
        for param in cur_params.keys():
            if param not in self._params:
                raise KeyError(f"Invalid param: {param}")
        self._params.update(cur_params)

    @property
    def data(self) -> list[T] | T:
        """Return the data of the current run."""
        return self._data

    @data.setter
    @abstractmethod
    def data(self) -> None:
        """Update the current data. Need to be further defined in subclass."""
        ...

    def precheck(self, force_rerun: bool = False) -> int:
        """
        Check the output folder and determine the rerun status.

        Rerun status:
        0: rerun all
        """
        rerun_status = 0
        if not self._output.exists():
            self._output.mkdir()

        if not self._output.is_dir():
            raise NotADirectoryError(f"{self._output} is already existed but not a folder. Aborted.")

        if force_rerun:
            shutil.rmtree(self._output)
            self._output.mkdir()

        if not any(self._output.iterdir()):
            pass
        elif self._ckp.exists():
            rerun_status = self._determine_rerun(*self.load_checkpoint())
        else:
            raise FileExistsError(f"Checkpoint file not found but the output directory {self._output} is not empty. Aborted.")

        return rerun_status

    def load_checkpoint(self) -> tuple[dict, ...]:
        """
        Load the checkpoint and retrieve the required params/data to determine the rerun status.

        This should be run before precheck.
        """
        logging.info("Loading from checkpoint...")
        try:
            with open(self._ckp, "rb") as f:
                params, *data = pickle.load(f)
        except (EOFError, ValueError):
            raise SystemExit("Align module checkpoint file corrupted. Please remove the output folder and try again.")

        return params, *data

    def save_checkpoint(self) -> None:
        """Save the parameters and data as a checkpoint for rerun."""
        if not hasattr(self, "_params"):
            raise AttributeError("No params available for checkpoint saving.")
        if not hasattr(self, "_data"):
            raise AttributeError("No data available for checkpoint saving.")
        else:
            if type(self._data) == list:
                if not all(self._data):
                    raise AttributeError("One item of the data is None.")
                data = self._data
            else:
                data = [self._data]
        with open(self._ckp, "wb") as f:
            pickle.dump((self._params, *data), f)

    @abstractmethod
    def _determine_rerun(self, prev_params: dict) -> int:
        """
        An abstract method for the actions that need to do when found a checkpoint file in the given output folder.

        Need to be further defined in the subclass.
        """
        logging.info("Get previous parameters and data and determine rerun status.")
        if not hasattr(self, "_params"):
            raise AttributeError(
                f"Currernt params not found. Please specify it through {self.__name__}.params = {dict} object first."
            )

    @staticmethod
    def _remove_files(files: Iterable[Path] | None = None, dirs: Iterable[Path] | None = None) -> None:
        """A static method to remove files and folders that need to be rerun."""
        if dirs:
            [shutil.rmtree(dir) for dir in dirs]
        if files:
            [file.unlink(missing_ok=True) for file in files]
