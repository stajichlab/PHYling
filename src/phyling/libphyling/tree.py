"""Filter and tree modules library."""

from __future__ import annotations

import gzip
import subprocess
import tempfile
from functools import wraps
from io import StringIO
from itertools import product
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Any, Callable, Literal, Sequence, TypeVar, overload

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

from .. import logger
from ..exception import BinaryNotFoundError
from ..external._libphykit import Saturation, compute_toverr
from . import SeqTypes, TreeMethods, TreeOutputFiles, _abc
from ._utils import CheckAttrs, CheckBinary, Timer, guess_seqtype, is_gzip_file

__all__ = ["MFA2Tree", "MFA2TreeList", "run_fasttree", "run_raxml", "run_iqtree"]
_C = TypeVar("Callable", bound=Callable[..., Any])


def _check_attributes(*attrs: str):
    """Decorator to ensure specific attributes are initialized before executing the function.

    Args:
        *attrs: Attribute names to check in the instance.

    Raises:
        AttributeError: If any specified attribute is `None` in the instance.
    """
    var_mapping = {
        "_tree": "build",
        "_toverr": "compute_toverr",
        "_saturation": "compute_saturation",
    }
    invalid_attrs = [attr for attr in attrs if attr not in var_mapping]
    if invalid_attrs:
        raise AttributeError(f"Invalid attribute names: {invalid_attrs}")

    def decorator(func: _C):
        @wraps(func)
        def wrapper(instance, *args, **kwargs):
            """Validate variable inequality and execute the wrapped function."""
            none_attrs = CheckAttrs.is_none(instance, *attrs)
            for var in sorted(none_attrs, key=lambda x: list(var_mapping.keys()).index(x)):
                raise AttributeError(f"Please run the {var_mapping[var]} method first.")
            return func(instance, *args, **kwargs)

        return wrapper

    return decorator


def _check_single_file(func: _C):
    """Decorator to ensure that the method is only called on a list with multiple MFA files.

    This decorator checks if the `MFA2TreeList` instance contains more than one MFA file. If there is only one file, it raises a
    `RuntimeError`. Otherwise, it calls the wrapped function.

    Args:
        func (Callable): The function to be wrapped.

    Returns:
        Callable: The decorated function.

    Raises:
        RuntimeError: If the `MFA2TreeList` instance contains fewer than two MFA files.
    """

    @wraps(func)
    def wrapper(self: MFA2TreeList, *args, **kwargs):
        if len(self) < 2:
            raise RuntimeError(f"{func.__name__} cannot be run since only found a single mfa.")
        return func(self, *args, **kwargs)

    return wrapper


class MFA2Tree(_abc.SeqFileWrapperABC):
    """Convert peptide/DNA multiple sequence alignment to Biopython Tree object.

    Optionally take a partition file to build tree with partition mode with raxml-ng and IQ-TREE.
    """

    __slots__ = ("_method", "_tree", "_toverr", "_saturation", "partition_file")

    @overload
    def __init__(self, file: str | Path) -> None: ...

    @overload
    def __init__(self, file: str | Path, name: str) -> None: ...

    def __init__(
        self,
        file: str | Path,
        name: str | None = None,
        *,
        partition_file: str | Path | None = None,
    ) -> None:
        """Initialize a MFA2Tree object.

        Initialize with a path of a multiple sequence alignment file (in fasta format) and a representative name (optional). The
        basename of the file will be used as the name if no name specified. The fasta won't be loaded immediately unless calling
        the load method.

        Args:
            file (str | Path): The path to the MSA file (plain or bgzipped).
            name (str, optional): A representative name for the sequences. Defaults to None.

        Examples:
            Load a peptide fasta:
            >>> MFA2Tree("data/101133at4751.aa.mfa")

            Load a peptide fasta with a given name:
            >>> MFA2Tree("data/101133at4751.aa.mfa", "msa_101133at4751")

            Load a bgzf compressed cds fasta:
            >>> MFA2Tree("data/101133at4751.cds.mfa", "msa_101133at4751")
        """
        super().__init__(file, name)
        self._method: str = None
        self._tree: Tree = None
        self._toverr: float = None
        self._saturation: float = None
        if partition_file:
            if isinstance(partition_file, (str, Path)):
                partition_file = Path(partition_file)
            else:
                raise TypeError(f"Argument partition_file only accepts str or Path. Got {type(partition_file)}.")
        self.partition_file: Path = partition_file

    @_abc.check_loaded
    def __len__(self) -> int:
        """Return the length of the MSA.

        Returns:
            int: The length of the MSA.
        """
        return self._data.get_alignment_length()

    @_abc.check_class
    def __gt__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is larger than another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is larger, otherwise False.
        """
        return self.toverr > other.toverr

    @_abc.check_class
    def __ge__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is larger than or equal to another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is larger or equal, otherwise False.
        """
        return self.toverr >= other.toverr

    @_abc.check_class
    def __lt__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is smaller than another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is smaller, otherwise False.
        """
        return self.toverr < other.toverr

    @_abc.check_class
    def __le__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is smaller than or equal to another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is smaller or equal, otherwise False.
        """
        return self.toverr <= other.toverr

    @_abc.check_class
    def __eq__(self, other: MFA2Tree) -> bool:
        """Check if the current object is equal to another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the objects have the same toverr, otherwise False.
        """
        return self.toverr == other.toverr

    def load(self) -> None:
        """Load the MultipleSeqAlignment.

        Assign the sequence type and its name after the MSA being loaded.
        """
        self._data = AlignIO.read(self.file, format="fasta")
        self._data.annotations["seqtype"] = self.seqtype
        self._data.annotations["seqname"] = self.name

    @property
    @_check_attributes("_tree")
    def tree(self) -> Tree:
        """Return the Biopython tree object.

        Raises:
            AttributeError: If tree haven't built.
        """
        return self._tree

    @property
    @_check_attributes("_tree")
    def method(self) -> str:
        """Return the method used for tree building.

        Raises:
            AttributeError: If tree haven't calculated.
        """
        return self._method

    @property
    @_check_attributes("_tree", "_toverr")
    def toverr(self) -> float:
        """Return the treeness/RCVs.

        Raises:
            AttributeError: If tree haven't built.
            AttributeError: If the toverr haven't calculated.
        """
        return self._toverr

    @property
    @_check_attributes("_tree", "_saturation")
    def saturation(self) -> float:
        """Return the saturation.

        Raises:
            AttributeError: If tree haven't built.
            AttributeError: If the saturation haven't calculated.
        """
        return self._saturation

    @overload
    def build(
        self,
        method: Literal["ft"],
        *,
        capture_cmd: bool = False,
    ) -> Tree: ...

    @overload
    def build(
        self,
        method: Literal["raxml"],
        output: str | Path | None = None,
        partition_file: str | Path | None = None,
        *,
        threads: int = 1,
        capture_cmd: bool = False,
    ) -> tuple[Tree, str]: ...

    @overload
    def build(
        self,
        method: Literal["iqtree"],
        output: str | Path | None = None,
        partition_file: str | Path | None = None,
        *,
        threads: int = 1,
        capture_cmd: bool = False,
    ) -> tuple[Tree, str]: ...

    @_abc.load_data
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        output: str | Path | None = None,
        partition_file: str | Path | None = None,
        *,
        threads: int = 1,
        capture_cmd: bool = False,
    ) -> Tree:
        """Build a phylogenetic tree using the specified method.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use raxml-ng.
                - "iqtree": Use IQ-TREE.
            output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
                Optional for raxml-ng and IQ-TREE.
            partition_file (str | Path | None, optional): Path to a partition file for model partitioning. Defaults to the
                object's partition_file if None.
            threads (int): Number of threads to use. Applicable for raxml-ng and IQ-TREE. Defaults to 1.
            capture_cmd (bool): If True, returns the command string used. Defaults to False.

        Returns:
            If capture_cmd is False, returns a Tree object. If capture_cmd is True, returns a tuple of the Tree object and the
            command string.
        """
        method = method.upper()
        if not partition_file and self.partition_file:
            partition_file = self.partition_file
        if partition_file and self.partition_file:
            logger.warning("MFA2Tree object has an already assigned partition_file. Use the one input to the build method.")
        cmd = None
        if method not in [m.name for m in TreeMethods]:
            raise KeyError(f'Invalid method: {method}. Expected one of: {", ".join([m.name.lower() for m in TreeMethods])}')
        if method == TreeMethods.FT.name:
            self._tree, cmd = run_fasttree(self, capture_cmd=True)
        elif method == TreeMethods.RAXML.name:
            self._tree, cmd = run_raxml(self, output, partition_file, threads=threads, capture_cmd=True)
        elif method == TreeMethods.IQTREE.name:
            self._tree, cmd = run_iqtree(self, output, partition_file, threads=threads, capture_cmd=True)
        else:
            raise NotImplementedError(f"Method not implemented yet: {TreeMethods[method.upper()].name}.")
        logger.debug("Tree building on %s is done.", self.name)
        self._method = method
        if capture_cmd:
            return self.tree, cmd
        else:
            return self.tree

    @_check_attributes("_tree")
    def compute_toverr(self) -> None:
        """Calculate the treeness/RCV of the tree by PhyKIT implementation."""
        self._toverr = compute_toverr(self.file, self.tree)

    @_check_attributes("_tree")
    def compute_saturation(self) -> None:
        """Calculate the saturation of the tree by PhyKIT implementation."""
        self._saturation = Saturation().compute_saturation(self._data, self.tree)

    def _guess_seqtype(self) -> str:
        """Guess and return the sequence type."""
        f = gzip.open(self.file, "rt") if is_gzip_file(self.file) else open(self.file)

        for r in FastaIO.SimpleFastaParser(f):
            seqtype = guess_seqtype(r[1], ignore_failed=True)
            if seqtype:
                f.close()
                return seqtype


class MFA2TreeList(_abc.SeqDataListABC[MFA2Tree]):
    """A wrapper for managing MFA2Tree objects, providing tree-building, concatenation and consensus tree calculations.

    This class allows working with multiple MFA2Tree instances, providing methods to build trees, compute metrics, and handle MSA
    concatenation efficiently with parallelization support.
    """

    _bound_class = MFA2Tree

    @overload
    def __init__(self, data: Sequence[str | Path | MFA2Tree]) -> None: ...

    @overload
    def __init__(self, data: Sequence[str | Path | MFA2Tree], names: Sequence[str]) -> None: ...

    def __init__(
        self,
        data: Sequence[str | Path | MFA2Tree] | None = None,
        names: Sequence[str] | None = None,
    ) -> None:
        """Initialize the MFA2TreeList wrapper object.

        Args:
            data (Sequence[str | Path | MFA2Tree] | None, optional): A sequence of data items.
            names (Sequence[str] | None, optional): A sequence of names corresponding to the data items.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to a MFA2Tree.
            KeyError: If the item already exists.
        """
        super().__init__(data, names)

    @overload
    def __getitem__(self, key: str) -> MFA2Tree: ...

    @overload
    def __getitem__(self, key: int) -> MFA2Tree: ...

    @overload
    def __getitem__(self, key: slice) -> MFA2TreeList: ...

    def __getitem__(self, key: str | int | slice) -> MFA2Tree | MFA2TreeList:
        """Retrieves an item or subset of items by name, index, or slice.

        Args:
            key (str | int | slice): The key to retrieve.

        Returns:
            MFA2Tree | MFA2TreeList: The corresponding item or subset of items.
        """
        return super().__getitem__(key)

    @property
    def trees(self) -> list[Tree]:
        """List of phylogenetic trees built from multiple sequence alignments.

        Returns:
            list[Tree]: A list of Biopython Tree objects.
        """
        return [mfa.tree for mfa in self]

    @property
    def method(self) -> str:
        """Return the method used for tree building.

        Raises:
            RuntimeError: If multiple methods are detected among the trees.

        Returns:
            str: The tree building method.
        """
        method = {mfa.method for mfa in self}
        if len(method) > 1:
            raise RuntimeError("Not all the trees built with the same method.")
        return next(iter(method))

    @property
    def toverrs(self) -> list[float]:
        """List of treeness/RCV values for each MFA2Tree.

        Returns:
            list[float]: Treeness/RCV values.
        """
        return [mfa.toverr for mfa in self]

    @property
    def saturations(self) -> list[float]:
        """List of saturation values for each MFA2Tree.

        Returns:
            list[float]: Saturation values.
        """
        return [mfa.saturation for mfa in self]

    def sort(self) -> None:
        """Sort the MFA2Tree objects by treeness/RCV in descending order."""
        super().sort(reverse=True)

    @Timer.timer
    def build(self, method: str, *, threads: int = 1, capture_cmd: bool = False, **kwargs) -> None:
        """Build trees for each MFA2Tree object.

        Args:
            method (str): Tree-building method to use.
            threads (int, optional): Number of parallel threads. Defaults to 1.
            capture_cmd (bool, optional): Capture command-line output. Defaults to False.
        """
        if len(self) == 1 or threads == 1:
            logger.debug("Sequential mode with maximum %s threads.", threads)
            for mfa2tree in self:
                _build_helper(
                    mfa2tree,
                    method,
                    kwargs.get("output"),
                    kwargs.get("partition_file"),
                    threads,
                    capture_cmd,
                )
        else:
            logger.debug("Multiprocesses mode: %s jobs are run concurrently.", threads)
            with ThreadPool(threads) as pool:
                pool.starmap(
                    _build_helper,
                    (
                        (
                            mfa2tree,
                            method,
                            kwargs.get("output"),
                            kwargs.get("partition_file"),
                            threads,
                            capture_cmd,
                        )
                        for mfa2tree in self
                    ),
                )

    @Timer.timer
    def compute_toverr(self, *, threads: int = 1) -> None:
        """Compute treeness/RCV for each tree.

        Args:
            threads (int, optional): Number of parallel threads. Defaults to 1.
        """
        logger.info("Calculating treeness/RCV...")
        if len(self) == 1 or threads == 1:
            logger.debug("Sequential mode with maximum %s threads.", threads)
            for mfa2tree in self:
                _compute_toverr_helper(mfa2tree)
        else:
            logger.debug("Multiprocesses mode: %s jobs are run concurrently.", threads)
            with ThreadPool(threads) as pool:
                pool.map(_compute_toverr_helper, self)
        logger.info("Done.")

    @Timer.timer
    def compute_saturation(self, *, threads: int = 1) -> None:
        """Compute saturation values for each tree.

        Args:
            threads (int, optional): Number of parallel threads. Defaults to 1.
        """
        logger.info("Calculating saturation...")
        if len(self) == 1 or threads == 1:
            logger.debug("Sequential mode with maximum %s threads.", threads)
            for mfa2tree in self:
                _compute_saturation_helper(mfa2tree)
        else:
            logger.debug("Multiprocesses mode: %s jobs are run concurrently.", threads)
            with ThreadPool(threads) as pool:
                pool.map(_compute_saturation_helper, self)
        logger.info("Done.")

    @Timer.timer
    @_check_single_file
    def get_consensus_tree(self) -> Tree:
        """Compute a consensus tree using ASTRAL.

        Returns:
            Tree: The consensus tree.

        Raises:
            BinaryNotFoundError: If ASTRAL is not installed.
            RuntimeError: If ASTRAL fails.
        """
        try:
            astral = CheckBinary.find("astral")
        except BinaryNotFoundError:
            raise BinaryNotFoundError(
                '%s not found. Please install it through "%s" or build the source following the instruction on %s'
                % (
                    "ASTRAL",
                    "conda install bioconda::aster",
                    "https://github.com/chaoszhang/ASTER",
                )
            )

        logger.info("Run ASTRAL to resolve consensus among multiple trees.")
        with tempfile.NamedTemporaryFile() as f:
            Phylo.write([mfa2tree.tree for mfa2tree in self], f.name, "newick")
            f.seek(0)
            cmd = [str(astral), f.name]
            try:
                result = subprocess.run(cmd, capture_output=True, check=True, text=True)
                logger.debug(
                    "Consensus tree inferred by %s is done:\n%s",
                    "ASTRAL",
                    result.stdout,
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"ASTRAL failed:\n{e.stdout}")
        return Phylo.read(StringIO(result.stdout), "newick")

    @overload
    def concat(self, output: Path, *, threads: int = 1) -> MFA2Tree: ...

    @overload
    def concat(self, output: Path, partition: bool, *, threads: int = 1) -> MFA2Tree: ...

    @Timer.timer
    @_check_single_file
    def concat(self, output: str | Path, partition: bool = False, *, threads: int = 1) -> MFA2Tree:
        """Concatenate selected MSAs into a single file and generate an MFA2Tree.

        Args:
            output (str | Path): Path to save the concatenated MSA file.
            partition (bool, optional): Whether to create a partition file. Defaults to False.
            threads (int, optional): Number of threads for parallel processing. Defaults to 1.

        Returns:
            MFA2Tree: An MFA2Tree object representing the concatenated result.

        Raises:
            TypeError: If the output is not a string or Path.
        """
        self.load()
        try:
            # Sort by MSA length with descending order (long -> short)
            if isinstance(output, (str, Path)):
                output = Path(output)
            else:
                raise TypeError("Argument output only accepts str or Path.")
            mfa2treelist = sorted(self, key=lambda msa: msa._data.get_alignment_length(), reverse=True)

            # Get sample names and create an empty MSA holder
            concat_alignments = MultipleSeqAlignment([])
            samples = set()
            msa_list = []
            for msa in mfa2treelist:
                msa_list.append(msa._data)
                for rec in msa._data:
                    samples.add(rec.id)
            samples = tuple(samples)
            for sample in samples:
                concat_alignments.append(SeqRecord(Seq(""), id=sample, description=""))
            concat_alignments.sort()
        finally:
            self.unload()

        # Fill missing taxon on some of the MSAs
        logger.info("Fill missing taxon...")
        if threads > 1:
            with ThreadPool(threads) as pool:
                msa_list = pool.starmap(_fill_missing_taxon, product([samples], msa_list))
        else:
            msa_list = [_fill_missing_taxon(samples, msa) for msa in msa_list]
        logger.info("Done.")

        # Concatenate
        logger.info("Concatenate selected MSAs...")

        start_idx = 0
        if partition:
            partition_info = []
        for msa in msa_list:
            # if self._get_missing_chars(mfa2tree):
            #     warnings.warn(f"{mfa2tree.name} contains missing state.", MissingStateWarning)
            #     continue
            concat_alignments += msa
            if partition:
                part_rec, end_idx = _generate_partition_record(msa, seqtype=self.seqtype, start=start_idx)
                partition_info.append(part_rec)
                start_idx = end_idx
        if partition:
            concat_alignments.annotations["partition"] = partition_info

        for seq in concat_alignments:
            seq.description = ""
        logger.info("Done.")

        output.mkdir(exist_ok=True)
        concat_file, partition_file = _output_concat_file(output, concat_alignments)
        return MFA2Tree(concat_file, partition_file=partition_file)


@_abc.check_loaded
def simple_distance_tree(mfa2tree: MFA2Tree, *, method: str) -> Tree:
    """[Deprecated] Run the tree calculation using a simple distance method."""
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, method)
    return constructor.build_tree(mfa2tree.load())


@overload
def run_fasttree(mfa2tree: MFA2Tree) -> Tree: ...


@overload
def run_fasttree(mfa2tree: MFA2Tree, *, capture_cmd: bool = True) -> tuple[Tree, str]: ...


@_abc.check_loaded
def run_fasttree(mfa2tree: MFA2Tree, *, capture_cmd: bool = False) -> Tree:
    """Runs FastTree to build a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree: The MFA2Tree object containing multiple sequence alignment data.
        capture_cmd: If True, returns the FastTree command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple containing the Tree object and the command string.
    """
    try:
        fasttree = CheckBinary.find(*TreeMethods.FT.bins)
    except BinaryNotFoundError:
        raise BinaryNotFoundError(
            '{} not found. Please install it through "{}"'.format(TreeMethods.FT.bins, "conda install bioconda::fasttree")
        )

    cmd = [str(fasttree), "-gamma", str(mfa2tree.file)]
    if mfa2tree.seqtype == SeqTypes.DNA:
        cmd.insert(1, "-gtr")
        cmd.insert(1, "-nt")
    else:
        cmd.insert(1, "-lg")

    try:
        logger.debug("Tree building cmd: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, check=True, text=True)
        logger.debug("Fasttree output:\n%s", result.stderr)  # This line output fasttree logs
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f'Tree building failed with cmd: {" ".join(cmd)}\n{e.stdout}')
    tree = Phylo.read(StringIO(result.stdout), "newick")
    if capture_cmd:
        return tree, " ".join(cmd)
    else:
        return tree


@overload
def run_raxml(
    mfa2tree: MFA2Tree,
    output: str | Path | None = None,
    partition_file: str | Path | None = None,
    *,
    threads: int = 1,
) -> Tree: ...


@overload
def run_raxml(
    mfa2tree: MFA2Tree,
    output: str | Path | None = None,
    partition_file: str | Path | None = None,
    *,
    threads: int = 1,
    capture_cmd: bool = True,
) -> tuple[Tree, str]: ...


@_abc.check_loaded
def run_raxml(
    mfa2tree: MFA2Tree,
    output: str | Path | None = None,
    partition_file: str | Path | None = None,
    *,
    threads: int = 1,
    capture_cmd: bool = False,
) -> Tree:
    """Runs raxml-ng to build a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree: The MFA2Tree object containing alignment data.
        output: Path to save the resulting tree file. If not provided, a temporary path is used.
        partition_file: Path to a partition file for model partitioning. Optional.
        threads: Number of threads to use for raxml-ng computation.
        capture_cmd: If True, returns the raxml-ng command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple of the Tree object and the command string.
    """
    try:
        raxml = CheckBinary.find(*TreeMethods.RAXML.bins)
    except BinaryNotFoundError:
        raise BinaryNotFoundError(
            '{} not found. Please install it through "{}"'.format(TreeMethods.RAXML.method, "conda install bioconda::raxml-ng")
        )

    if output:
        if not isinstance(output, (str, Path)):
            raise TypeError(f"Argument output only accepts str or Path. Got {type(output)}.")
    if partition_file:
        if not isinstance(partition_file, (str, Path)):
            raise TypeError(f"Argument partition only accepts str or Path. Got {type(output)}.")
        partition_file = Path(partition_file)
        if not partition_file.exists():
            raise FileNotFoundError(f"{partition_file}")
        if not partition_file.is_file():
            raise RuntimeError(f"{partition_file} is not a file.")

    tempdir = None
    try:
        if not output:
            tempdir = tempfile.TemporaryDirectory()
            output = tempdir.name

        if not isinstance(output, (str, Path)):
            raise TypeError(f"Argument output only accepts str or Path. Got {type(output)}.")
        output = Path(output) / TreeMethods.RAXML.method.lower()
        output.mkdir(exist_ok=True)

        cmd = [
            str(raxml),
            "--all",
            "--msa",
            str(mfa2tree.file),
            "--prefix",
            str(output / mfa2tree.name),
            "--bs-trees",
            "1000",
            "--threads",
            str(threads),
        ]
        if partition_file:
            cmd.extend(["--model", str(partition_file)])
        elif mfa2tree.seqtype == SeqTypes.DNA:
            cmd.extend(["--model", "GTR+G"])
        else:
            cmd.extend(["--model", "LG+G"])

        try:
            logger.debug("Tree building cmd: %s", " ".join(cmd))
            result = subprocess.run(cmd, capture_output=True, check=True, text=True)
            logger.debug("%s output:\n%s", TreeMethods.RAXML.method, result.stdout)  # This line output raxml logs
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f'Tree building failed with cmd: {" ".join(cmd)}\n{e.stdout}')
        tree = Phylo.read(output / f"{mfa2tree.name}.raxml.support", "newick")
        if capture_cmd:
            return tree, " ".join(cmd)
        else:
            return tree

    finally:
        if tempdir:
            tempdir.cleanup()


@overload
def run_iqtree(
    mfa2tree: MFA2Tree,
    output: str | Path | None = None,
    partition_file: str | Path | None = None,
    *,
    threads: int = 1,
) -> Tree: ...


@overload
def run_iqtree(
    mfa2tree: MFA2Tree,
    output: str | Path | None = None,
    partition_file: str | Path | None = None,
    *,
    threads: int = 1,
    capture_cmd: bool = True,
) -> tuple[Tree, str]: ...


@_abc.check_loaded
def run_iqtree(
    mfa2tree: MFA2Tree,
    output: str | Path | None = None,
    partition_file: str | Path | None = None,
    *,
    threads: int = 1,
    capture_cmd: bool = False,
) -> Tree | tuple[Tree, str]:
    """Runs IQ-TREE to construct a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree: The MFA2Tree object containing alignment data.
        output: Path to save the resulting tree file. If not provided, a temporary path is used.
        partition_file: Path to a partition file for model partitioning. Optional.
        threads: Number of threads to use for IQ-TREE computation.
        capture_cmd: If True, returns the IQ-TREE command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple of the Tree object and the command string.
    """
    try:
        iqtree = CheckBinary.find(*TreeMethods.IQTREE.bins)
    except BinaryNotFoundError:
        raise BinaryNotFoundError(
            '{} not found. Please install it through "{}"'.format(TreeMethods.IQTREE.method, "conda install bioconda::iqtree")
        )

    if output:
        if not isinstance(output, (str, Path)):
            raise TypeError(f"Argument output only accepts str or Path. Got {type(output)}.")
    if partition_file:
        if not isinstance(partition_file, (str, Path)):
            raise TypeError(f"Argument partition only accepts str or Path. Got {type(output)}.")
        partition_file = Path(partition_file)
        if not partition_file.exists():
            raise FileNotFoundError(f"{partition_file}")
        if not partition_file.is_file():
            raise RuntimeError(f"{partition_file} is not a file.")

    tempdir = None
    try:
        if not output:
            tempdir = tempfile.TemporaryDirectory()
            output = tempdir.name

        output = Path(output) / TreeMethods.IQTREE.method.lower()
        output.mkdir(exist_ok=True)

        cmd = [
            str(iqtree),
            "-s",
            str(mfa2tree.file),
            "--prefix",
            str(output / mfa2tree.name),
            "-bb",
            "1000",
            "-T",
            "AUTO",
            "-ntmax",
            str(threads),
        ]
        if mfa2tree.seqtype == SeqTypes.DNA:
            cmd.extend(["-m", "GTR+G"])
        else:
            cmd.extend(["-m", "LG+G"])
        if partition_file:
            if not isinstance(partition_file, (str, Path)):
                raise TypeError(f"Argument partition only accepts str or Path. Got {type(output)}.")
            cmd.extend(["-p", str(partition_file)])

        try:
            logger.debug("Tree building cmd: %s", " ".join(cmd))
            result = subprocess.run(cmd, capture_output=True, check=True, text=True)
            logger.debug("%s output:\n%s", TreeMethods.IQTREE.method, result.stdout)  # This line output iqtree logs
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f'Tree building failed with cmd: {" ".join(cmd)}\n{e.stdout}')
        tree = Phylo.read(output / f"{mfa2tree.name}.contree", "newick")
        if capture_cmd:
            return tree, " ".join(cmd)
        else:
            return tree

    finally:
        if tempdir:
            tempdir.cleanup()


def _compute_toverr_helper(instance: MFA2Tree) -> None:
    """Helper function to run the `compute_toverr` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `compute_toverr` method.
    """
    instance.compute_toverr()


def _compute_saturation_helper(instance: MFA2Tree) -> None:
    """Helper function to run the `compute_saturation` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `compute_saturation` method.
    """
    instance.compute_saturation()


def _build_helper(
    instance: MFA2Tree,
    method: Literal["ft", "raxml", "iqtree"],
    output: str | Path | None,
    partition_file: str | Path | None,
    threads: int,
    verbose: bool,
) -> Tree:
    """Helper function to run the `build` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `build` method.
        method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use raxml-ng.
                - "iqtree": Use IQ-TREE.
        output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
                Optional for raxml-ng and IQ-TREE.
        partition_file (str | Path | None, optional): Path to a partition file for model partitioning. Defaults to the
            object's partition_file if None.
        threads (int): Number of threads to use. Applicable for raxml-ng and IQ-TREE. Defaults to 1.
        capture_cmd (bool): If True, returns the command string used. Defaults to False.

    Returns:
        Tree: The resulting phylogenetic tree.
    """
    instance.build(
        method,
        output=output,
        partition_file=partition_file,
        threads=threads,
        capture_cmd=verbose,
    )


def _fill_missing_taxon(samples: Sequence[str], msa: MultipleSeqAlignment) -> MultipleSeqAlignment:
    """Add gap-only sequences to the alignment for missing taxa.

    Args:
        samples (Sequence[str]): List of all expected taxa.
        msa (MultipleSeqAlignment): The existing multiple sequence alignment.

    Returns:
        MultipleSeqAlignment: The updated alignment with added gap-only sequences for missing taxa.
    """
    missing = set(samples) - {seq.id for seq in msa}
    for sample in missing:
        msa.append(
            SeqRecord(
                Seq("-" * msa.get_alignment_length()),
                id=sample,
                description="",
            )
        )
        msa.sort()
    return msa


def _generate_partition_record(msa: MultipleSeqAlignment, *, seqtype, start: int = 0) -> tuple[str, int]:
    """Generate a partition record for raxml-ng or IQ-TREE.

    Args:
        msa (MultipleSeqAlignment): The multiple sequence alignment.
        seqtype: The type of sequence (DNA or protein).
        start (int, optional): The starting position for the partition. Defaults to 0.

    Returns:
        tuple[str, int]: A tuple containing the partition record string and the end position.
    """

    # Assign model for raxml/iqtree compatible partition file
    model = "GTR+G" if seqtype == SeqTypes.DNA else "LG+G"
    end = msa.get_alignment_length()
    return (
        f'{model}, {msa.annotations["seqname"]} = {start + 1}-{start + end}',
        start + end,
    )


def _output_concat_file(output: Path, concat_alignments: MultipleSeqAlignment) -> tuple[Path, Path]:
    """Write concatenated FASTA and partition files to the specified output directory.

    Args:
        output (Path): The directory where files will be written.
        concat_alignments (MultipleSeqAlignment): The concatenated multiple sequence alignment.

    Returns:
        tuple[Path, Path | None]: Paths to the concatenated FASTA file and partition file, or None if partitioning is disabled.
    """
    concat_file = output / TreeOutputFiles.CONCAT
    with open(concat_file, "w") as f:
        SeqIO.write(concat_alignments, f, format="fasta")
        logger.info("Concatenated fasta is output to %s", concat_file)

    if "partition" in concat_alignments.annotations:
        partition_file = output / TreeOutputFiles.PARTITION
        with open(partition_file, "w") as f:
            for part_rec in concat_alignments.annotations["partition"]:
                f.write(f"{part_rec}\n")
        logger.info("Partition file is output to %s", partition_file)
    else:
        partition_file = None
        logger.debug("Partition mode disabled. No partition file output.")
    return concat_file, partition_file
