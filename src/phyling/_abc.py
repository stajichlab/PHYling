import logging
import pickle
import shutil
import textwrap
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, Iterable, Iterator, Sequence, TypeVar

import phyling.config as config
import phyling.exception as exception
from phyling._utils import check_cls_vars

T = TypeVar("T")
SFW = TypeVar("SFW", bound="SeqFileWrapperABC")


class SeqFileWrapperABC(ABC):
    """An abstract class for sequence files."""

    class Decorator:
        """A decorator to contains classmethods "check_seqtypes"."""

        @classmethod
        def check_class(cls, func: callable):
            """Check whether self and other are the same class."""

            def wrapper(self: SFW, other: SFW):
                if not isinstance(other, self.__class__):
                    raise TypeError(f"Can only operate with another {self.__class__} object but got {type(other)}")
                return func(self, other)

            return wrapper

        @classmethod
        def check_seqtype(cls, func: callable):
            """Check whether self and other represent the same seqtype."""

            def wrapper(self: SFW, other: SFW):
                if self.seqtype != other.seqtype:
                    raise exception.SeqtypeError("Items represent different seqtypes.")
                return func(self, other)

            return wrapper

    def __init__(self, file: str | Path, *, seqtype: str | None = None) -> None:
        """Initialize the object and check sequence type."""
        self._file = Path(file)
        if not self._file.exists():
            raise FileNotFoundError(f"{self._file}")
        if self._file.is_dir():
            raise IsADirectoryError(f"{self._file}")
        if seqtype:
            self._set_seqtype(seqtype)

    def __repr__(self) -> str:
        """Return the string representation."""
        return f"{self.__class__.__qualname__}({self.name}" + ", %s)"

    @abstractmethod
    @Decorator.check_seqtype
    @Decorator.check_class
    def __eq__(self, other: SFW) -> bool:
        """Return true if the two objects have the same value."""
        ...

    @abstractmethod
    @Decorator.check_seqtype
    @Decorator.check_class
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
    def __init__(self, files: Iterator[str | Path | T] | Sequence[str | Path | T]) -> None:
        """Initialize the object and store data into a list."""
        self.data: list[T] = []
        if not isinstance(files, Iterator | tuple | list):
            raise TypeError(f"{self.__class__.__qualname__} is not an iterator/tuple/list")

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


class OutputPrecheckABC(ABC):
    """A abstract class that provides features for input/output precheck and checkpoint loading/saving."""

    folder: Path
    ckp: str

    def __new__(cls, *args, **kwargs):
        """Prevent instantiation by raising an exception."""
        raise NotImplementedError("This class cannot be instantiated.")

    @classmethod
    @abstractmethod
    def setup(cls, folder: Path = None, ckp: str = None) -> None:
        """Setup the class variable."""
        if folder is not None:
            cls.folder = folder
        if ckp is not None:
            cls.ckp = ckp
        for var in ("folder", "ckp"):
            if getattr(cls, var, None) is None:
                raise ValueError(f"The class variable '{var}' must be set before calling this method.")

    @classmethod
    @check_cls_vars("folder", "ckp")
    @abstractmethod
    def precheck(cls, params: dict, data: tuple, force_rerun: bool = False) -> tuple | None:
        """
        Check the output folder and determine the rerun status.

        Rerun status:
        0: rerun all
        """
        if not cls.folder.exists():
            cls.folder.mkdir()

        if not cls.folder.is_dir():
            raise NotADirectoryError(f"{cls.folder} is already existed but not a folder. Aborted.")

        if force_rerun:
            shutil.rmtree(cls.folder)
            cls.folder.mkdir()

        if any(cls.folder.iterdir()) and not (cls.folder / cls.ckp).exists():
            raise FileExistsError(f"Checkpoint file not found but the output directory {cls.folder} is not empty. Aborted.")

    @classmethod
    @check_cls_vars("folder", "ckp")
    def load_checkpoint(cls) -> tuple[dict, ...]:
        """
        Load the checkpoint and retrieve the required params/data to determine the rerun status.

        This should be run before output precheck.
        """
        logging.info("Loading from checkpoint...")
        try:
            with open(cls.folder / cls.ckp, "rb") as f:
                params, *data = pickle.load(f)
        except FileNotFoundError:
            raise FileNotFoundError("Checkpoint file not found.")
        except (EOFError, ValueError):
            raise RuntimeError("Align module checkpoint file corrupted. Please remove the output folder and try again.")

        return params, *data

    @classmethod
    @check_cls_vars("folder", "ckp")
    def save_checkpoint(cls, params: dict, *data) -> None:
        """Save the parameters and data as a checkpoint for rerun."""
        with open(cls.folder / cls.ckp, "wb") as f:
            pickle.dump((params, *data), f)

    @classmethod
    @check_cls_vars("folder", "ckp")
    @abstractmethod
    def _determine_rerun(cls, cur: tuple, prev: tuple) -> tuple:
        """
        An abstract method for the actions that need to do when found a checkpoint file in the given output folder.

        Need to be further defined in the subclass.
        """
        ...

    @staticmethod
    def _remove_files(files: Iterable[Path] | None = None, dirs: Iterable[Path] | None = None) -> None:
        """A static method to remove files and folders that need to be rerun."""
        if dirs:
            [shutil.rmtree(dir) for dir in dirs]
        if files:
            [file.unlink(missing_ok=True) for file in files]
