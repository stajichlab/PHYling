"""Filter and tree modules library."""

from __future__ import annotations

import gzip
import tempfile
import traceback
from functools import wraps
from itertools import product
from multiprocessing import Manager
from multiprocessing.pool import ThreadPool
from multiprocessing.sharedctypes import Synchronized
from multiprocessing.synchronize import Condition
from pathlib import Path
from typing import Callable, Literal, Sequence, TypeVar, overload

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

from .. import logger
from ..external import (
    Astral,
    Concordance,
    FastTree,
    Iqtree,
    ModelFinder,
    PartitionRecord,
    Partitions,
    Raxml,
    RaxmlHandler,
    UFBoot,
)
from ..external._libphykit import Saturation, compute_toverr
from . import SeqTypes, TreeMethods, TreeOutputFiles, _abc
from ._utils import CheckAttrs, Timer, guess_seqtype, is_gzip_file, progress_daemon

__all__ = ["MFA2Tree", "MFA2TreeList"]
_R = TypeVar("R")


def _check_attributes(*attrs: str) -> Callable[[Callable[..., _R]], Callable[..., _R]]:
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

    def decorator(func: Callable[..., _R]) -> Callable[..., _R]:
        @wraps(func)
        def wrapper(instance, *args, **kwargs):
            """Validate variable inequality and execute the wrapped function."""
            none_attrs = CheckAttrs.is_none(instance, *attrs)
            for var in sorted(none_attrs, key=lambda x: list(var_mapping.keys()).index(x)):
                raise AttributeError(f"Please run the {var_mapping[var]} method first.")
            return func(instance, *args, **kwargs)

        return wrapper

    return decorator


def _check_single_file(func: Callable[..., _R]) -> Callable[..., _R]:
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

    __slots__ = ("_method", "_tree", "_cmds", "_toverr", "_saturation")

    @overload
    def __init__(self, file: str | Path, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...

    @overload
    def __init__(self, file: str | Path, name: str, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...

    def __init__(self, file: str | Path, name: str | None = None, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None:
        """Initialize a MFA2Tree object.

        Initialize with a path of a multiple sequence alignment file (in fasta format) and a representative name (optional). The
        basename of the file will be used as the name if no name specified. The fasta won't be loaded immediately unless calling
        the load method.

        Args:
            file (str | Path): The path to the MSA file (plain or bgzipped).
            name (str, optional): A representative name for the sequences. Defaults to None.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Examples:
            Load a peptide fasta:
            >>> MFA2Tree("data/101133at4751.aa.mfa")

            Load a peptide fasta with a given name:
            >>> MFA2Tree("data/101133at4751.aa.mfa", "msa_101133at4751")

            Load a bgzf compressed cds fasta:
            >>> MFA2Tree("data/101133at4751.cds.mfa", "msa_101133at4751")
        """
        super().__init__(file, name, seqtype=seqtype)
        self._method: str = None
        self._tree: Tree = None
        self._cmds: list[str] = []
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
    def cmds(self) -> list[str]:
        """Return the commands that were used in the build method.

        Raises:
            AttributeError: If tree haven't built.
        """
        return self._cmds

    @property
    @_check_attributes("_tree")
    def method(self) -> str:
        """Return the method used for tree building.

        Raises:
            AttributeError: If tree haven't built.
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
        output: str | Path | None = None,
        model: str | Path = "AUTO",
        *,
        noml: bool = False,
        seed: int = -1,
    ): ...

    @overload
    def build(
        self,
        method: Literal["raxml", "iqtree"],
        output: str | Path | None = None,
        model: str | Path = "AUTO",
        *,
        bs: int = 0,
        scfl: int = 0,
        seed: int = -1,
        threads: int = 1,
        threads_max: int = 1,
    ): ...

    @_abc.load_data
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        output: str | Path | None = None,
        model: str | Path = "AUTO",
        *,
        noml: bool = False,
        bs: int = 0,
        scfl: int = 0,
        seed: int = -1,
        threads: int = 1,
        threads_max: int = 1,
        **kwargs,
    ) -> Tree:
        """Build a phylogenetic tree using the specified method.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
            output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
                Optional for RAxML-NG and IQ-TREE.
            model (str | Path | None, optional): String of model or path to a partition file for model partitioning. Defaults to
                "AUTO" that use ModelFinder to identify the best model.
            noml (bool): Disable Maximum-likelihood estimation to speed up tree building. Only available for FastTree. Defaults to
                False.
            bs (int): UFBoot value. Defaults to 0 (disabled).
            scfl (int): Site concordance factor threshold for IQ-TREE. Must be at least 1000 if set. Defaults to 0 (disabled).
            threads (int): Number of threads to use. Applicable for RAxML-NG and IQ-TREE. Defaults to 1.

        Returns:
            Tree: A Biopython Tree object.
        """
        method = method.upper()
        if method not in [m.name for m in TreeMethods]:
            raise KeyError(f"Invalid method: {method}. Expected one of: {', '.join([m.name.lower() for m in TreeMethods])}")

        tempdir = None
        try:
            self._cmds.clear()
            if not output:
                tempdir = tempfile.TemporaryDirectory()
                output = tempdir.name

            output = Path(output)
            output.mkdir(exist_ok=True)

            if model == "AUTO":
                modelfinder = ModelFinder(
                    self.file,
                    output / "modelfinder",
                    seqtype=self.seqtype,
                    method=method.lower(),
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                )
                logger.info("Find the best-fit model...")
                modelfinder.run()
                model = modelfinder.result
                self._cmds.append(f"ModelFinder cmd: {modelfinder.cmd}")
                logger.info(f"Best-fit model: {model}")
                logger.info("Done.")
            elif Path(model).is_file():
                # Partitioning analysis
                if method == TreeMethods.FT.name:
                    raise ValueError(f"{TreeMethods.FT.method} does not support partitioning analysis.")

            if method == TreeMethods.FT.name:
                runner = FastTree(
                    self.file,
                    output / f"{self.file.name}.nw",
                    seqtype=self.seqtype,
                    model=model,
                    noml=noml,
                    **kwargs,
                )
            elif method == TreeMethods.RAXML.name:
                runner = Raxml(
                    self.file,
                    output / self.file.name,
                    seqtype=self.seqtype,
                    model=model,
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                    **kwargs,
                )
            elif method == TreeMethods.IQTREE.name:
                runner = Iqtree(
                    self.file,
                    output / self.file.name,
                    seqtype=self.seqtype,
                    model=model,
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                    **kwargs,
                )
            else:
                raise NotImplementedError(f"Method not implemented yet: {TreeMethods[method].name}.")

            logger.info(f"Build tree by {runner._prog}...")
            runner.run()
            tree_file = runner.result
            self._method = method
            model = runner.model
            self._cmds.append(f"Tree building cmd: {runner.cmd}")
            logger.info("Done.")

            if bs > 0:
                runner = UFBoot(self.file, tree_file, output / "ufboot", model=model, bs=bs, seed=seed, threads=threads)
                logger.info("Bootstrapping...")
                runner.run()
                tree_file = runner.result
                self._cmds.append(f"UFBoot cmd: {runner.cmd}")
                logger.info("Done.")

            if scfl > 0:
                runner = Concordance(self.file, tree_file, output / "concord", model=model, scfl=scfl, seed=seed, threads=threads)
                logger.info("Calculate site concordance factor...")
                runner.run()
                tree_file = runner.result
                self._cmds.append(f"Branch concordance cmd: {runner.cmd}")
                logger.info("Done.")

            self._tree = Phylo.read(tree_file, "newick")
            logger.debug("Tree building on %s is done.", self.name)

        finally:
            if tempdir:
                tempdir.cleanup()

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
    __slots__ = ("_cmds",)

    @overload
    def __init__(self, data: Sequence[str | Path | MFA2Tree], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...

    @overload
    def __init__(
        self, data: Sequence[str | Path | MFA2Tree], names: Sequence[str], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO"
    ) -> None: ...

    def __init__(
        self,
        data: Sequence[str | Path | MFA2Tree] | None = None,
        names: Sequence[str] | None = None,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ) -> None:
        """Initialize the MFA2TreeList wrapper object.

        Args:
            data (Sequence[str | Path | MFA2Tree] | None, optional): A sequence of data items.
            names (Sequence[str] | None, optional): A sequence of names corresponding to the data items.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to a MFA2Tree.
            KeyError: If the item already exists.
        """
        self._cmds: list[str] = []
        super().__init__(data, names, seqtype=seqtype)

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
    def cmds(self) -> list[str]:
        """Return the commands that were used in the build method."""
        return [cmd for mfa2tree in self for cmd in mfa2tree.cmds] + self._cmds

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
        noml: bool = False,
        bs: int = 0,
        scfl: int = 0,
        threads: int = 1,
    ) -> None:
        """Build trees for each MFA2Tree object.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
            noml (bool): Disable Maximum-likelihood estimation to speed up tree building. Only available for FastTree. Defaults to
                False.
            bs (int): UFBoot value. Defaults to 0 (disabled).
            scfl (int): Site concordance factor threshold for IQ-TREE. Must be at least 1000 if set. Defaults to 0 (disabled).
            threads (int, optional): Number of parallel threads. Defaults to 1.
        """
        with Manager() as manager:
            counter = manager.Value("i", 0)
            condition = manager.Condition()

            progress = progress_daemon(len(self), counter, condition, step=min(max(10, len(self) // 200 * 50), 500))
            progress.start()

            params_per_task = [(mfa2tree, method, None, model, noml, bs, scfl, threads, counter, condition) for mfa2tree in self]
            try:
                if len(self) == 1 or threads == 1:
                    logger.debug("Sequential mode with %s threads.", threads)
                    for params in params_per_task:
                        _build_helper(*params)
                else:
                    logger.debug("Multiprocesses mode with %s jobs and 1 thread for each.", threads)
                    with ThreadPool(threads) as pool:
                        pool.starmap(_build_helper, params_per_task)
            except Exception:
                logger.error("%s", traceback.format_exc())
                raise

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
    def get_consensus_tree(self, *, seed: int = -1, threads: int = 1) -> Tree:
        """Compute a consensus tree using ASTRAL.

        Returns:
            Tree: The consensus tree.

        Raises:
            BinaryNotFoundError: If ASTRAL is not installed.
            RuntimeError: If ASTRAL fails.
        """
        with tempfile.NamedTemporaryFile() as f_in, tempfile.NamedTemporaryFile() as f_out:
            self._cmds.clear()
            Phylo.write([mfa2tree.tree for mfa2tree in self], f_in.name, "newick")
            f_in.seek(0)
            runner = Astral(f_in.name, f_out.name, seed=seed, threads=threads)
            runner.run()
            self._cmds.append(runner.cmd)
            return Phylo.read(runner.result, "newick")

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
            mfa2treelist = sorted(self, key=lambda msa: len(msa), reverse=True)

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
        partition_info = Partitions("raxml")
        model = "GTR" if self.seqtype == SeqTypes.DNA else "LG"
        for msa in msa_list:
            concat_alignments += msa
            end_idx = start_idx + msa.get_alignment_length()
            partition_info.add(
                PartitionRecord(
                    name=msa.annotations["seqname"],
                    start=start_idx + 1,
                    end=end_idx,
                    model=model,
                )
            )
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
        with RaxmlHandler(partition_file, "w") as f:
            f.write(partition_info)
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
    noml: bool = False,
    bs: int = 0,
    scfl: int = 0,
    threads: int = 1,
    counter: Synchronized | None = None,
    condition: Condition | None = None,
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
        noml (bool): Disable Maximum-likelihood estimation to speed up tree building. Only available for FastTree. Defaults to
            False.
        bs (int): UFBoot value. Defaults to 0 (disabled).
        scfl (int): Site concordance factor threshold for IQ-TREE. Must be at least 1000 if set. Defaults to 0 (disabled).
        threads (int): Number of threads to use. Applicable for RAxML-NG and IQ-TREE. Defaults to 1.

    Returns:
        Tree: The resulting phylogenetic tree.
    """
    instance.build(method, output, model, noml=noml, bs=bs, scfl=scfl, threads=threads)

    if counter and condition:
        with condition:
            counter.value += 1
            condition.notify()


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
