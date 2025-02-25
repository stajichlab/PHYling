"""Filter and tree modules library."""

from __future__ import annotations

import gzip
import tempfile
from functools import wraps
from itertools import product
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Any, Callable, Literal, Sequence, TypeVar, overload

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

from .. import logger
from ..external import Astral, Concordance, FastTree, Iqtree, ModelFinder, Raxml, UFBoot
from ..external._libphykit import Saturation, compute_toverr
from . import SeqTypes, TreeMethods, TreeOutputFiles, _abc
from ._utils import CheckAttrs, Timer, guess_seqtype, is_gzip_file

__all__ = ["MFA2Tree", "MFA2TreeList"]
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


def _check_single_file(func: _C) -> _C:
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

    Optionally take a partition file to build tree with partition mode with RAxML-NG and IQ-TREE.
    """

    __slots__ = ("_method", "_tree", "_toverr", "_saturation")

    @overload
    def __init__(self, file: str | Path) -> None: ...

    @overload
    def __init__(self, file: str | Path, name: str) -> None: ...

    def __init__(self, file: str | Path, name: str | None = None) -> None:
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
        method: Literal["ft", "raxml", "iqtree"],
        output: str | Path,
        model: str | Path = "AUTO",
        *,
        bs: int = 1000,
        scfl: int = 100,
        threads: int = 1,
        capture_cmd: bool = True,
    ) -> tuple[Tree, str]: ...

    @_abc.load_data
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        output: str | Path | None = None,
        model: str | Path = "AUTO",
        *,
        bs: int = 1000,
        scfl: int = 100,
        threads: int = 1,
        capture_cmd: bool = False,
    ) -> Tree:
        """Build a phylogenetic tree using the specified method.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
            output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
                Optional for RAxML-NG and IQ-TREE.
            model (str | Path | None, optional): Path to a partition file for model partitioning. Defaults to the
                object's partition_file if None.
            bs (int): Bootstrap value. Defaults to 100.
            threads (int): Number of threads to use. Applicable for RAxML-NG and IQ-TREE. Defaults to 1.
            capture_cmd (bool): If True, returns the command string used. Defaults to False.

        Returns:
            If capture_cmd is False, returns a Tree object. If capture_cmd is True, returns a tuple of the Tree object and the
            command string.
        """
        method = method.upper()
        tempdir = None
        capture_cmd_msg = ""
        try:
            if not output:
                tempdir = tempfile.TemporaryDirectory()
                output = tempdir.name

            output = Path(output)
            output.mkdir(exist_ok=True)

            if model == "AUTO":
                modelfinder = ModelFinder(
                    self.file, output / "modelfinder", seqtype=self.seqtype, method=method.lower(), threads=threads
                )
                model = modelfinder()
            elif Path(model).is_file():
                # Partitioning analysis
                if method == TreeMethods.FT.name:
                    raise ValueError(f"{TreeMethods.FT.method} does not support partitioning analysis.")

            if method not in [m.name for m in TreeMethods]:
                raise KeyError(f"Invalid method: {method}. Expected one of: {', '.join([m.name.lower() for m in TreeMethods])}")
            if method == TreeMethods.FT.name:
                runner = FastTree(self.file, output / f"{self.file.name}.nw", seqtype=self.seqtype, model=model, threads=threads)
            elif method == TreeMethods.RAXML.name:
                runner = Raxml(self.file, output / self.file.name, seqtype=self.seqtype, model=model, threads=threads)
            elif method == TreeMethods.IQTREE.name:
                runner = Iqtree(self.file, output / self.file.name, seqtype=self.seqtype, model=model, threads=threads)
            else:
                raise NotImplementedError(f"Method not implemented yet: {TreeMethods[method].name}.")
            tree_file = runner()
            self._method = method
            model = runner.model
            capture_cmd_msg += f"Tree building cmd: {runner.cmd}\n"

            if bs > 0:
                runner = UFBoot(self.file, tree_file, output / "ufboot", model=model, bs=bs, threads=threads)
                tree_file = runner()
                capture_cmd_msg += f"UFBoot cmd: {runner.cmd}\n"

            if scfl > 0:
                runner = Concordance(self.file, tree_file, output / "concord", model=model, scfl=scfl, threads=threads)
                tree_file = runner()
                capture_cmd_msg += f"Concordance factor cmd: {runner.cmd}\n"

            self._tree = Phylo.read(tree_file, "newick")
            logger.debug("Tree building on %s is done.", self.name)

        finally:
            if tempdir:
                tempdir.cleanup()

        if capture_cmd:
            return self.tree, capture_cmd_msg
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
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        model: str | Path = "AUTO",
        *,
        bs: int = 0,
        scfl: int = 0,
        threads: int = 1,
        capture_cmd: bool = False,
    ) -> None:
        """Build trees for each MFA2Tree object.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
            bs (int): Bootstrap value. Defaults to 100.
            threads (int, optional): Number of parallel threads. Defaults to 1.
            capture_cmd (bool, optional): Capture command-line output. Defaults to False.
        """
        if len(self) == 1 or threads == 1:
            logger.debug("Sequential mode with maximum %s threads.", threads)
            for mfa2tree in self:
                _build_helper(
                    mfa2tree,
                    method,
                    None,
                    model,
                    bs,
                    scfl,
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
                            None,
                            model,
                            bs,
                            scfl,
                            1,
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
    def get_consensus_tree(self, *, threads: int = 1) -> Tree:
        """Compute a consensus tree using ASTRAL.

        Returns:
            Tree: The consensus tree.

        Raises:
            BinaryNotFoundError: If ASTRAL is not installed.
            RuntimeError: If ASTRAL fails.
        """
        with tempfile.NamedTemporaryFile() as f_in, tempfile.NamedTemporaryFile() as f_out:
            Phylo.write([mfa2tree.tree for mfa2tree in self], f_in.name, "newick")
            f_in.seek(0)
            runner = Astral(f_in.name, f_out.name, threads=threads)
            return Phylo.read(runner(), "newick")

    @Timer.timer
    @_check_single_file
    def concat(self, output: str | Path, *, threads: int = 1) -> tuple[Path, Path]:
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
        partition_info = []
        for msa in msa_list:
            concat_alignments += msa
            part_rec, end_idx = _generate_partition_record(msa, seqtype=self.seqtype, start=start_idx)
            partition_info.append(part_rec)
            start_idx = end_idx

        for seq in concat_alignments:
            seq.description = ""
        logger.info("Done.")

        output.mkdir(exist_ok=True)
        concat_file = output / TreeOutputFiles.CONCAT
        with open(concat_file, "w") as f:
            SeqIO.write(concat_alignments, f, format="fasta")
            logger.info("Concatenated fasta is output to %s", concat_file)

        partition_file = output / TreeOutputFiles.PARTITION
        with open(partition_file, "w") as f:
            for part_rec in partition_info:
                f.write(f"{part_rec}\n")
        logger.info("Partition file is output to %s", partition_file)
        return concat_file, partition_file


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
    output: str | Path | None = None,
    model: str | Path = "AUTO",
    bs: int = 1000,
    scfl: int = 100,
    threads: int = 1,
    capture_cmd: bool = False,
) -> Tree:
    """Helper function to run the `build` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `build` method.
        method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
        output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
                Optional for RAxML-NG and IQ-TREE.
        partition_file (str | Path | None, optional): Path to a partition file for model partitioning. Defaults to the
            object's partition_file if None.
        bs (int): Bootstrap value. Defaults to 100.
        threads (int): Number of threads to use. Applicable for RAxML-NG and IQ-TREE. Defaults to 1.
        capture_cmd (bool): If True, returns the command string used. Defaults to False.

    Returns:
        Tree: The resulting phylogenetic tree.
    """
    instance.build(method, output, model, bs=bs, scfl=scfl, threads=threads, capture_cmd=capture_cmd)


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
    """Generate a partition record for RAxML-NG or IQ-TREE.

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
        f"{model}, {msa.annotations['seqname']} = {start + 1}-{start + end}",
        start + end,
    )
