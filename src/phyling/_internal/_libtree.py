"""Phylogenetic tree construction methods."""

from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import hashlib
import re
import shutil
import subprocess
import tempfile
from io import StringIO
from itertools import product
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Callable, Iterator, Sequence

import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import NCBICodonTable, NCBICodonTableDNA
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import phyling._internal._abc as _abc
import phyling._internal._config as _config
import phyling._internal._utils as _utils
import phyling.exception as exception
from phyling import logger
from phyling._internal._libphykit import Saturation, compute_toverr
from phyling._internal._libphyling import OutputPrecheck as libphylingPrecheck

expected_files = {
    "msas_dir": "selected_MSAs",
    "treeness": "treeness.tsv",
    "concat": f"concat_alignments.{_config.aln_ext}",
    "partition": "concat_alignments.partition",
    "tree_nw": "final_tree.nw",
    "tree_img": "final_tree.png",
    "raxml": "raxml",
    "iqtree": "iqtree",
}


class MFA2Tree(_abc._SeqFileWrapperABC):
    """
    Convert peptide/DNA multiple sequence alignment to Biopython Tree object and calculate the treeness/RCV.

    Optionally take a partition file to build tree with partition mode with RAxML and IQTree.
    """

    class Decorator:
        """A decorator to contains classmethods "check_built" and "check_toverr"."""

        @classmethod
        def check_built(cls, func: Callable):
            """Check whether the method build() has already being run."""

            def wrapper(self: MFA2Tree, *args, **kwargs):
                if not hasattr(self, "_tree"):
                    raise AttributeError("build() need to be run first.")
                return func(self, *args, **kwargs)

            return wrapper

        @classmethod
        def check_toverr(cls, func: Callable):
            """Check whether the method compute_toverr() has already being run."""

            def wrapper(self: MFA2Tree, *args, **kwargs):
                if not hasattr(self, "_toverr"):
                    raise AttributeError("compute_toverr() need to be run first.")
                return func(self, *args, **kwargs)

            return wrapper

    def __init__(self, file: str | Path, partition: Path | None = None, *, seqtype: str | None = None) -> None:
        """Initialize the object and check sequence type."""
        super().__init__(file, seqtype=seqtype)
        if partition:
            self._partition = partition
        if not seqtype:
            self._seqtype = None
            self._determine_seqtype()

    def __repr__(self) -> str:
        """Return the string representation."""
        info = f"seqtype={self.seqtype}"
        if hasattr(self, "_toverr"):
            info += f", toverr={self.toverr}"
        return super().__repr__() % info

    def __eq__(self, other: MFA2Tree) -> bool:
        """Return true if the two objects have the same toverr."""
        super().__eq__(other)
        return self.toverr == other.toverr

    def __lt__(self, other: MFA2Tree) -> bool:
        """Return true if the toverr of the current object is smaller than another object."""
        super().__lt__(other)
        return self.toverr < other.toverr

    @property
    def name(self) -> str:
        """Return the basename of the file."""
        return self.path.name

    @property
    def checksum(self) -> int:
        """Compute the checksum of the stored data."""
        return _utils.get_file_checksum(self.path)

    @property
    def msa(self) -> MultipleSeqAlignment:
        """Return Biopython MultipleSeqAlignment object."""
        alignment = AlignIO.read(self.path, format="fasta")
        alignment.annotations["seqtype"] = self._seqtype
        alignment.annotations["seqname"] = self.name
        return alignment

    @property
    @Decorator.check_built
    def tree(self) -> Tree:
        """Return the Biopython tree object."""
        return self._tree

    @property
    @Decorator.check_built
    def method(self) -> str:
        """Return the treeness/RCVs."""
        return self._method

    @property
    @Decorator.check_built
    @Decorator.check_toverr
    def toverr(self) -> float:
        """Return the treeness/RCVs."""
        return self._toverr

    def build(self, method: str, output: Path | None = None, *, threads: int = 1, verbose=False) -> None:
        """Run tree building step and return the Biopython Tree object."""
        if method not in _config.avail_tree_methods.keys():
            raise KeyError(f'Argument "method" not falling in {_config.avail_tree_methods.keys()}')
        if method in ("upgma", "nj"):
            tree = self._with_phylo_module(method=method)
        elif method == "ft":
            tree = self._with_FastTree(verbose=verbose)
        elif method == "raxml":
            tree = self._with_RAxML(output=output, threads=threads, verbose=verbose)
        elif method == "iqtree":
            tree = self._with_IQTree(output=output, threads=threads, verbose=verbose)
        else:
            raise NotImplementedError(f"{_config.avail_tree_methods[method]} is not implemented yet.")
        logger.debug(f"Tree building on {self.name} is done.")
        self._method = method
        self._tree = tree

    @Decorator.check_built
    def compute_toverr(self) -> None:
        """Calculate the treeness/RCV of the tree by PhyKIT implementation."""
        self._toverr = compute_toverr(self.path, self.tree)

    @Decorator.check_built
    def compute_saturation(self) -> None:
        """Calculate the saturation of the tree by PhyKIT implementation."""
        self._saturation = Saturation.compute_saturation()

    def _determine_seqtype(self) -> str:
        chars = set()
        [chars.update(str(rec.seq)) for rec in self.msa]
        if "-" in chars:
            chars.remove("-")
        if len(chars) <= 4 and chars.issubset(set(NCBICodonTableDNA.nucleotide_alphabet)):
            self._seqtype = _config.seqtype_cds
        elif chars.issubset(set(NCBICodonTable.protein_alphabet)):
            self._seqtype = _config.seqtype_pep
        else:
            raise exception.SeqtypeError("Cannot determine seqtype. Aborted.")

    def _with_phylo_module(self, *, method: str) -> Tree:
        """[Deprecated] Run the tree calculation using a simple distance method."""
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, method)
        return constructor.build_tree(self.msa)

    def _with_FastTree(self, *, verbose: bool) -> Tree:
        """Run the tree calculation using FastTree."""
        fasttree_bins = ("FastTree", "fasttree")

        for fasttree in fasttree_bins:
            fasttree = shutil.which(fasttree)
            if fasttree:
                if verbose:
                    logger.debug(f"Found FastTree at {fasttree}")
                break

        if not fasttree:
            raise exception.BinaryNotFoundError(
                'FastTree not found. Please install it through "conda install bioconda::fasttree"'
            )

        cmd = [fasttree, "-gamma", str(self.path)]
        if self.seqtype == _config.seqtype_cds:
            cmd.insert(1, "-nt")
        else:
            cmd.insert(1, "-lg")
        if verbose:
            logger.info(f'Tree building cmd: {" ".join(cmd)}')

        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        tree, err = p.communicate()
        if not tree:
            raise RuntimeError(err)
        return Phylo.read(StringIO(tree), "newick")

    def _with_RAxML(self, output: Path, *, threads: int, verbose: bool) -> Tree:
        """Run the tree calculation using RAxML."""
        raxml_bins = (
            "raxmlHPC-PTHREADS-AVX2",
            "raxmlHPC-PTHREADS-SSE3",
            "raxmlHPC-PTHREADS",
            "raxmlHPC-AVX2",
            "raxmlHPC-SSE3",
            "raxmlHPC",
        )
        for raxml in raxml_bins:
            raxml = shutil.which(raxml)
            if raxml:
                if verbose:
                    logger.debug(f"Found RAxML at {raxml}")
                break

        if not raxml:
            raise exception.BinaryNotFoundError('RAxML not found. Please install it through "conda install bioconda::raxml"')

        if output is None:
            tempdir = Path(tempfile.mkdtemp())
            output_raxml = tempdir / "raxml"
        else:
            output_raxml = output / "raxml"

        output_raxml.mkdir(exist_ok=True)
        cmd = [
            raxml,
            "-s",
            str(self.path),
            "-p",
            "12345",
            "-w",
            str((output_raxml).absolute()),
            "-n",
            self.seqtype,
        ]
        if self.seqtype == _config.seqtype_cds:
            cmd.extend(["-m", "GTRCAT"])
        else:
            cmd.extend(["-m", "PROTGAMMAAUTO"])
        if hasattr(self, "_partition"):
            cmd.extend(["-q", str(self._partition)])
            with open(self._partition) as f:
                line_counts = f.read().count("\n")
            threads = max(min(line_counts, threads), 2)
        if re.match(r".*PTHREADS.*", raxml):
            cmd.extend(["-T", str(threads)])
        if verbose:
            logger.info(f'Tree building cmd: {" ".join(cmd)}')

        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stdout, _ = p.communicate()
        try:
            tree = Phylo.read(output_raxml / f"RAxML_bestTree.{self.seqtype}", "newick")
        except FileNotFoundError:
            raise RuntimeError(f"RAxML run on {str(self.path)} failed: {stdout}")
        if output is None:
            shutil.rmtree(tempdir)

        return tree

    def _with_IQTree(self, output: Path, *, threads: int, verbose: bool) -> Tree:
        """Run the tree calculation using IQTree."""
        iqtree_bins = ("iqtree", "iqtree2")

        for iqtree in iqtree_bins:
            iqtree = shutil.which(iqtree)
            if iqtree:
                if verbose:
                    logger.debug(f"Found IQTree at {iqtree}")
                break

        if not iqtree:
            raise exception.BinaryNotFoundError('IQTree not found. Please install it through "conda install bioconda::iqtree"')

        if output is None:
            tempdir = Path(tempfile.mkdtemp())
            output_iqtree = tempdir / "iqtree"
        else:
            output_iqtree = output / "iqtree"

        output_iqtree.mkdir(exist_ok=True)
        cmd = [
            iqtree,
            "-s",
            str(self.path),
            "--prefix",
            str(output_iqtree / "iqtree"),
            "-m",
            "MFP",
            "--mset",
            "raxml",
            "-T",
            "AUTO",
            "-ntmax",
            str(threads),
        ]
        if hasattr(self, "_partition"):
            cmd.extend(["-p", str(self._partition)])
        if verbose:
            logger.info(f'Tree building cmd: {" ".join(cmd)}')

        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stdout, stderr = p.communicate()
        if stderr:
            raise RuntimeError(f"IQTree run on {str(self.path)} failed: {stderr}")

        try:
            tree = Phylo.read(output_iqtree / "iqtree.treefile", "newick")
        except FileNotFoundError:
            raise RuntimeError(f"IQTree run on {str(self.path)} success but not able to load tree from file: {stdout}")

        if output is None:
            shutil.rmtree(tempdir)

        return tree


class MFA2TreeWrapper(_abc._DataListABC[MFA2Tree]):
    """A Wrapper for MFA2Tree object(s) which also provides concatenated data preparation and consensus tree calculation."""

    class Decorator:
        """A decorator to contains classmethods "check_single_file" and "check_sorted"."""

        @classmethod
        def check_single_file(cls, func: Callable):
            """Check whether data has more than one MFA2Tree object and raise an error if not."""

            def wrapper(self: MFA2TreeWrapper, *args, **kwargs):
                if len(self) < 2:
                    raise RuntimeError(f"{func.__name__} cannot be run since only single mfa found after filtering by toverr.")
                return func(self, *args, **kwargs)

            return wrapper

        @classmethod
        def check_sorted(cls, func: Callable):
            """Check whether data is sorted and raise an error if not."""

            def wrapper(self: MFA2TreeWrapper, *args, **kwargs):
                results = func(self, *args, **kwargs)
                if self.is_sorted:
                    return results
                if isinstance(results, MFA2Tree):
                    return results
                raise RuntimeError("Data is not sorted. Please run compute_toverr() to sort it by toverr first.")

            return wrapper

    def __init__(self, files: list[str | Path | MFA2Tree], *, seqtype: str | None = None) -> None:
        """Initialize the wrapper object."""
        super().__init__(files)
        for file in files:
            if isinstance(file, (str, Path)):
                data = MFA2Tree(file, seqtype=seqtype)
            elif isinstance(file, MFA2Tree):
                data = file
            else:
                raise TypeError(f"{self.__class__.__qualname__} only accepts list of str/Path/MFA2Tree")
            if not hasattr(self, "_seqtype"):
                self._seqtype = data.seqtype
            self._append(data)
        self._issorted = False

    def __repr__(self) -> str:
        """Return the string representation."""
        info = f"seqtype={self.seqtype}, sorted={self.is_sorted}"
        if self.is_sorted:
            info += f", method={self.method}"
        return super().__repr__() % f"({info})"

    @Decorator.check_sorted
    def __getitem__(self, idx: int | slice) -> MFA2Tree | MFA2TreeWrapper[MFA2Tree]:
        """Get the data by index or slice. Note that the data need to be sorted before getting with slice."""
        data = self.data[idx]
        if isinstance(data, MFA2Tree):
            return data
        elif isinstance(data, list):
            wrapper = self.__class__(data, seqtype=self.seqtype)
            wrapper.sort()
            return wrapper
        else:
            raise TypeError("Index must be an integer or slice")

    @property
    def seqtype(self) -> str:
        """Return whether the sequences are peptide or DNA."""
        return self._seqtype

    @property
    def checksum(self) -> str:
        """Compute md5 checksum based on the crc checksum of individual files."""
        checksums = 0
        for data in self:
            checksums += data.checksum
        return hashlib.md5(str(checksums).encode()).hexdigest()

    @property
    def method(self) -> str:
        """Return the method that was used for tree building."""
        method = {mfa.method for mfa in self}
        if len(method) > 1:
            raise RuntimeError("Not all the trees built with the same method.")
        return next(iter(method))

    @property
    def is_sorted(self) -> bool:
        """Return whether the data is sorted."""
        return self._issorted

    def sort(self) -> None:
        """Sort the MFA2Tree objects by treeness/RCV."""
        super().sort(reverse=True)
        self._issorted = True

    @_utils.timing
    def build(self, method: str, output: Path | None = None, *, threads: int = 1, verbose=False) -> None:
        """Run the tree building step for each MFA2Tree object."""
        if len(self) == 1 or threads == 1:
            logger.debug(f"Sequential mode with maximum {threads} threads.")
            for mfa2tree in self:
                self._build_helper(mfa2tree, method=method, output=output, threads=threads, verbose=verbose)
        else:
            logger.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
            with ThreadPool(threads) as pool:
                pool.starmap(
                    self._build_helper,
                    ((mfa2tree, method, output, threads, verbose) for mfa2tree in self),
                )

    @_utils.timing
    def compute_toverr(self, *, threads: int = 1) -> None:
        """Calculate the treeness/RCV for each MFA2Tree object and sort the MFA2Tree object by the values."""
        logger.info("Calculating treeness/RCV...")
        if len(self) == 1 or threads == 1:
            logger.debug(f"Sequential mode with maximum {threads} threads.")
            for mfa2tree in self:
                self._compute_toverr_helper(mfa2tree)
        else:
            logger.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
            with ThreadPool(threads) as pool:
                pool.map(self._compute_toverr_helper, self)
        self.sort()
        logger.info("Done.")

    @_utils
    def compute_saturation(self, *, threads: int = 1) -> None:
        """Calculate the saturation for each MFA2Tree object."""
        logger.info("Calculating saturation...")
        if len(self) == 1 or threads == 1:
            logger.debug(f"Sequential mode with maximum {threads} threads.")
            for mfa2tree in self:
                self._compute_saturation_helper(mfa2tree)
        else:
            logger.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
            with ThreadPool(threads) as pool:
                pool.map(self._compute_saturation_helper, self)
        logger.info("Done.")

    @_utils.timing
    @Decorator.check_sorted
    @Decorator.check_single_file
    def get_consensus_tree(self) -> Tree:
        """Run astral to get the consensus tree."""
        if not shutil.which("astral"):
            raise exception.BinaryNotFoundError(
                "ASTRAL not found. "
                'Please install it through "conda install bioconda::aster" '
                "or build the source following the instruction on https://github.com/chaoszhang/ASTER"
            )

        logger.info("Run ASTRAL to resolve consensus among multiple trees.")
        temp = StringIO()
        Phylo.write([mfa2tree.tree for mfa2tree in self], temp, "newick")
        temp.seek(0)
        p = subprocess.Popen(
            ["astral", "/dev/stdin"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, _ = p.communicate(temp.read())
        temp.close()
        logger.info("Done.")
        return Phylo.read(StringIO(stdout), "newick")

    @_utils.timing
    @Decorator.check_sorted
    @Decorator.check_single_file
    def concat(
        self, output: Path, samples: Iterator[str] | Sequence[str] | None, *, partition: str | None = None, threads: int = 1
    ) -> MFA2Tree:
        """Concatenate selected MSAs to a single mfa file and return the corresponding MFA2Tree object."""
        alignmentList = [mfa2tree.msa for mfa2tree in self]
        # Sort by MSA length, since shorter MSA might have higher chances to have missing states
        alignmentList = sorted(alignmentList, key=lambda alignment: alignment.get_alignment_length())

        # Create an empty MSA holder
        concat_alignments = MultipleSeqAlignment([])
        if not samples:
            samples = set()
            for mfa2tree in self:
                for rec in mfa2tree.msa:
                    samples.add(rec.id)
        for sample in samples:
            concat_alignments.append(SeqRecord(Seq(""), id=sample, description=""))
        concat_alignments.sort()

        # Fill missing taxon on some of the MSAs
        logger.info("Fill missing taxon...")
        if threads > 1:
            with ThreadPool(threads) as pool:
                alignmentList = pool.starmap(self._fill_missing_taxon, product([samples], alignmentList))
        else:
            alignmentList = [self._fill_missing_taxon(samples, alignment) for alignment in alignmentList]
        logger.info("Done.")

        # Save partition_info
        if partition:
            logger.info("Check msa with missing states...")
            alignmentList = self._concat_msa_with_missing_states(alignmentList)
            logger.info("Done.")
            if len(alignmentList) == 1:
                raise RuntimeError(
                    "Only 1 MSA left after combining msas with missing states. Please rerun with --partition disabled."
                )

        # Concatenate
        logger.info("Concatenate selected MSAs...")
        for alignment in alignmentList:
            concat_alignments += alignment
        if partition:
            concat_alignments.annotations["partition"] = self._get_partition_info(alignmentList, partition)
        for seq in concat_alignments:
            seq.description = ""
        logger.info("Done.")

        concat_file, partition_file = self._output_concat_file(output, concat_alignments)

        return MFA2Tree(concat_file, partition=partition_file, seqtype=self.seqtype)

    @staticmethod
    def _fill_missing_taxon(samples: Iterator[str] | Sequence[str], alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
        """Include empty gap-only strings to the alignment for taxa lacking an ortholog."""
        missing = set(samples) - {seq.id for seq in alignment}
        for sample in missing:
            alignment.append(
                SeqRecord(
                    Seq("-" * alignment.get_alignment_length()),
                    id=sample,
                    description="",
                )
            )
            alignment.sort()
        return alignment

    def _append(self, item: MFA2Tree) -> None:
        """Add new sample to the collection."""
        if not isinstance(item, MFA2Tree):
            raise TypeError(f"Can only append {MFA2Tree.__qualname__} object but got {type(item)}.")
        if item.seqtype != self.seqtype:
            raise exception.SeqtypeError("Item represents different seqtype.")
        if item.name in self._names():
            raise exception.IdenticalKeyError(f"File {item.name} already exists.")
        self.data.append(item)

    def _build_helper(self, mfa2tree: MFA2Tree, method: str, output: Path | None, threads: int, verbose: bool) -> Tree:
        """Helper function to run the MFA2tree build method."""
        mfa2tree.build(method, output, threads=threads, verbose=verbose)

    def _compute_toverr_helper(self, mfa2tree: MFA2Tree) -> None:
        """Helper function to run the MFA2Tree compute_toverr method."""
        mfa2tree.compute_toverr()

    def _compute_saturation_helper(self, mfa2tree: MFA2Tree) -> None:
        """Helper function to run the MFA2Tree compute_toverr method."""
        mfa2tree.compute_saturation()

    def _concat_msa_with_missing_states(self, alignmentList: list[MultipleSeqAlignment]) -> list[MultipleSeqAlignment]:
        """Concatenate the msas with missing states to prevent warning from RAxML."""
        idx = 0
        count = 0
        alignmentList_length_orig = len(alignmentList)
        makeup_partition = MultipleSeqAlignment([], annotations={"seqtype": self.seqtype})
        for rec in alignmentList[0]:
            makeup_partition.append(SeqRecord(Seq(""), id=rec.id, description=""))
        while idx < len(alignmentList):
            alignment = alignmentList[idx]
            missing_chars = self._get_missing_chars(alignment)
            if missing_chars:
                makeup_partition += alignment
                alignmentList.pop(idx)
                count += 1
            else:
                idx += 1

        if self._get_missing_chars(makeup_partition) and alignmentList:
            makeup_partition += alignmentList.pop(0)
            count += 1

        makeup_partition.annotations["seqname"] = f"Concatenated_partition_with_{count}_MSAs"
        msg = f"{count} / {alignmentList_length_orig} of MSAs have missing states."
        if count / alignmentList_length_orig >= 0.5:
            logger.warning(msg)
            logger.warning(
                "This indicates samples might contain many short sequences. Recommended to run with --partition disabled."
            )
        else:
            logger.info(msg)
        logger.info("MSAs with missing states would be concatenated into a single partition.")

        alignmentList.insert(0, makeup_partition)

        return alignmentList

    def _get_missing_chars(self, alignment: MultipleSeqAlignment) -> set | None:
        """Get the missing states by lookup the DNA or peptide codon table."""
        chars = set()
        for rec in alignment:
            chars.update(set(rec.seq))

        if "-" in chars:
            chars.remove("-")

        All_chars = (
            set(NCBICodonTable.protein_alphabet)
            if self.seqtype == _config.seqtype_pep
            else set(NCBICodonTableDNA.nucleotide_alphabet)
        )
        missing_chars = All_chars - chars

        return missing_chars if missing_chars else None

    def _get_partition_info(self, alignmentList: list[MultipleSeqAlignment], partition: str):
        """Get raxml-supported partition info."""
        end = 0
        partition_info = []

        if partition in ("codon", "seq+codon") and self.seqtype == _config.seqtype_pep:
            logger.warning('Peptide sequence detected. Force to use "seq" mode for partitioning.')
            OutputPrecheck.partition = partition = "seq"
        if partition:
            logger.info(f"Partitioning with {partition} mode.")

        # Assign model for raxml/iqtree compatible partition file
        model = "DNA" if self.seqtype == _config.seqtype_cds else "AUTO"
        for alignment in alignmentList:
            start = end + 1
            end += alignment.get_alignment_length()
            if partition == "seq":
                partition_info.append(f'{model}, {alignment.annotations["seqname"]} = {start}-{end}')
            elif partition == "seq+codon":
                for pos in range(3):
                    partition_info.append(f'{model}, {alignment.annotations["seqname"]}_codon{pos + 1} = {start + pos}-{end}\\3')
        if partition == "codon":
            for pos in range(3):
                partition_info.append(f"{model}, codon{pos + 1} = {pos + 1}-{end}\\3")
        return partition_info

    def _output_concat_file(self, output: Path, concat_alignments: MultipleSeqAlignment) -> tuple[Path, Path]:
        """Output the concatenated fasta and partition file."""
        concat_file = output / expected_files["concat"]
        with open(concat_file, "w") as f:
            SeqIO.write(concat_alignments, f, format="fasta")
            logger.info(f"Concatenated fasta is output to {concat_file}")
        if "partition" in concat_alignments.annotations:
            partition_file = output / expected_files["partition"]
            with open(partition_file, "w") as f:
                [print(line, file=f) for line in concat_alignments.annotations["partition"]]
            logger.info(f"Partition file is output to {partition_file}")
        else:
            partition_file = None
            logger.debug("The concat_alignments object doesn't have partition info. No partition file output.")
        return concat_file, partition_file


class OutputPrecheck(_abc._OutputPrecheckABC):
    """A class that provides features for input/output precheck, checkpoint loading/saving and final tree output."""

    folder: Path
    ckp: str
    method: str
    concat: bool
    partition: str
    figure: bool

    @classmethod
    def setup(cls, *, folder: Path = None, ckp: str = _config.phylotree_checkpoint, **kwargs) -> None:
        """Setup the class variable."""
        super().setup(folder, ckp)
        for key, value in kwargs.items():
            if key in ("method", "concat", "partition", "figure"):
                setattr(cls, key, value)

    @classmethod
    @_utils.check_cls_vars("folder", "ckp")
    def precheck(cls, params: dict, *, force_rerun: bool = False) -> tuple[int, MFA2TreeWrapper | None]:
        """
        Check the output folder and determine the rerun status.

        0: rerun all
        1: rerun treeness/RCV filtering
        2: rerun concatenate step
        3: rerun final tree building step
        """
        if params.keys() != _config.phylotree_precheck_params:
            raise KeyError(f"Params should contain keys {_config.phylotree_precheck_params}")
        super().precheck(params, (), force_rerun=force_rerun)
        if not any(cls.folder.iterdir()):
            return (0, None)
        results = cls._determine_rerun(params, *cls.load_checkpoint())
        return results

    @classmethod
    def load_checkpoint(cls) -> tuple[dict, MFA2TreeWrapper]:
        """
        Load the checkpoint and retrieve the required params/data to determine the rerun status.

        This should be run before precheck.
        """
        return super().load_checkpoint()

    @classmethod
    def save_checkpoint(cls, params: dict, wrapper: MFA2TreeWrapper) -> None:
        """Save the parameters and data as a checkpoint for rerun."""
        params[f'{expected_files["msas_dir"]}_checksum'] = _utils.get_multifiles_checksum(
            (cls.folder / expected_files["msas_dir"]).glob(f"*.{_config.aln_ext}")
        )
        super().save_checkpoint(params, wrapper)

    @classmethod
    @_utils.check_cls_vars("folder", "method", "concat", "partition", "figure")
    def output_results(cls, tree: Tree):
        """Output the tree in newick format and figure."""
        Phylo.draw_ascii(tree)
        output_tree = cls.folder / expected_files["tree_nw"]
        logger.info(f"Output tree to {output_tree}")
        with open(output_tree, "w") as f:
            strategy = "concatenate" if cls.concat else "consensus"
            print(
                f"# Final tree is built using {_config.avail_tree_methods[cls.method]} with {strategy} strategy",
                file=f,
            )
            if cls.concat and cls.partition:
                print(f"# Partition is enabled using {cls.partition} mode", file=f)
            Phylo.write(tree, f, "newick")

        if cls.figure:
            fig, ax = plt.subplots(figsize=(20, 12))
            output_fig = cls.folder / expected_files["tree_img"]
            logger.info(f"Output figure to {output_fig}")
            Phylo.draw(tree, axes=ax)
            fig.savefig(output_fig)

    @classmethod
    def _determine_rerun(cls, cur_params: dict, prev_params: dict, wrapper: MFA2TreeWrapper) -> tuple[int, MFA2TreeWrapper]:
        """Define the actions that need to do when found a checkpoint file in the given output folder."""

        def rerun_status_3(rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
            """
            Return rerun status as 3 and the corresponding folders and files need to be removed.

            rm_dirs: all folders related to raxml and iqtree
            rm_files: all final tree files (nw and png, if any)
            """
            for x in cls.folder.iterdir():
                if x.is_dir():
                    if x.name != expected_files["msas_dir"]:
                        rm_dirs.add(x)
                else:
                    if x.name in (expected_files["tree_nw"], expected_files["tree_img"]):
                        rm_files.add(cls.folder / x)
            return 3, rm_dirs, rm_files

        def rerun_status_2(rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
            """
            Return rerun status as 2 and the corresponding folders and files need to be removed.

            rm_dirs: none
            rm_files: all files except checkpoint and treeness
            """
            _, rm_dirs, rm_files = rerun_status_3(rm_dirs, rm_files)
            rm_files.update(
                x
                for x in cls.folder.iterdir()
                if x.is_file() and x.name not in (expected_files["treeness"], _config.phylotree_checkpoint)
            )

            return 2, rm_dirs, rm_files

        def rerun_status_1(rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
            """
            Return rerun status as 1 and the corresponding folders and files need to be removed.

            rm_dirs: all folders
            rm_files: all files except checkpoint
            """
            _, rm_dirs, rm_files = rerun_status_2(rm_dirs, rm_files)
            for x in cls.folder.iterdir():
                if x.is_dir():
                    rm_dirs.add(x)
                else:
                    if x.name != _config.phylotree_checkpoint:
                        rm_files.add(x)
            return 1, rm_dirs, rm_files

        def rerun_status_0(rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
            """
            Return rerun status as 0 and the corresponding folders and files need to be removed.

            rm_dirs: all folders
            rm_files: all files
            """
            _, rm_dirs, rm_files = rerun_status_1(rm_dirs, rm_files)
            for x in cls.folder.iterdir():
                if x.is_dir():
                    rm_dirs.add(x)
                else:
                    rm_files.add(x)

            return 0, rm_dirs, rm_files

        mfas = [file for file in (cls.folder / expected_files["msas_dir"]).glob(f"*.{_config.aln_ext}")]
        if mfas:
            mfas_checksum = _utils.get_multifiles_checksum(mfas)
        else:
            mfas_checksum = -1
        cur_params.update({f'{expected_files["msas_dir"]}_checksum': mfas_checksum})
        diff_params = {param[0] for param in set(cur_params.items()) ^ set(prev_params.items())}
        if not diff_params:
            raise SystemExit("Files not changed and parameters are identical to the previous run. Aborted.")
        logger.debug(f"{cls._determine_rerun.__name__}.diff_params = {diff_params}")
        if "inputs" in diff_params:
            raise SystemExit(f"Inputs are changed but the output directory {cls.folder} is not empty. Aborted.")
        if not cur_params["concat"]:
            print("Rerun all. Remove all previous files.")
            rerun_status, rm_dirs, rm_files = rerun_status_0()
        else:
            if "concat" in diff_params:
                logger.info("Rerun all. Remove all previous files.")
                rerun_status, rm_dirs, rm_files = rerun_status_0()
            elif diff_params.intersection((f'{expected_files["msas_dir"]}_checksum', "top_n_toverr")):
                logger.info("Rerun treeness/RCV filtering.")
                rerun_status, rm_dirs, rm_files = rerun_status_1()
            elif "partition" in diff_params:
                logger.info("Rerun concatenate step.")
                rerun_status, rm_dirs, rm_files = rerun_status_2()
            elif "method" in diff_params:
                logger.info("Rerun final tree building step.")
                rerun_status, rm_dirs, rm_files = rerun_status_3()

        cls._remove_files(files=rm_files, dirs=rm_dirs)
        if rerun_status == 0:
            wrapper = None
        return rerun_status, wrapper


def single_mfa(
    input: Path,
    output: Path,
    *,
    method: str,
    seqtype: str | None = None,
    threads: int = 1,
) -> tuple[Tree, MFA2TreeWrapper]:
    """
    Run the analysis in single mfa mode.

    Return the tree and wrapper for rerun purpose.
    """
    wrapper = MFA2TreeWrapper(input, seqtype=seqtype)
    wrapper.build(method=method, output=output, threads=threads, verbose=True)
    wrapper.compute_toverr(threads=threads)
    _output_treeness(wrapper, output)

    logger.info("Only one MSA is selected. Report the tree directly.")
    return wrapper[0].tree, wrapper


def concat(
    inputs: list,
    output: Path,
    *,
    method: str,
    seqtype: str | None = None,
    samples: list | None = None,
    top_n_toverr: int = 50,
    partition: str | None = None,
    threads: int = 1,
    rerun: int = 0,
    wrapper: MFA2TreeWrapper | None = None,
) -> tuple[Tree, MFA2TreeWrapper]:
    """
    Run the analysis in concat mode.

    Required more than 1 files. Return the tree and wrapper for rerun purpose.
    """
    if bool(rerun) ^ bool(wrapper):
        raise RuntimeError("Rerun must to be coupled with wrapper.")

    if rerun == 0:
        wrapper = MFA2TreeWrapper(inputs, seqtype=seqtype)
        logger.info("Use FastTree to generate trees for treeness filtering.")
        wrapper.build(method="ft", threads=threads)
        wrapper.compute_toverr(threads=threads)
    else:
        wrapper = wrapper

    if rerun <= 1:
        _output_treeness(wrapper, output, n=top_n_toverr)

    if 0 < top_n_toverr < len(wrapper):
        corrected_top_n = len(wrapper[:top_n_toverr])
        logger.info(f"Use the MSAs which the treeness are among the top {corrected_top_n} for final tree building.")
        wrapper = wrapper[:corrected_top_n]
    else:
        logger.info("Use all MSAs for final tree building.")

    if rerun <= 2:
        concat_mfa2tree = wrapper.concat(output=output, samples=samples, partition=partition, threads=threads)
    else:
        partition_file = output / expected_files["partition"] if partition else None
        concat_mfa2tree = MFA2Tree(output / expected_files["concat"], partition=partition_file, seqtype=seqtype)

    logger.info("Use the concatednated fasta to generate final tree.")
    concat_mfa2tree.build(method=method, output=output, threads=threads, verbose=True)
    return concat_mfa2tree.tree, wrapper


def consensus(
    inputs: list,
    output: Path,
    *,
    method: str,
    seqtype: str | None = None,
    top_n_toverr: int = 50,
    threads: int = 1,
) -> tuple[Tree, MFA2TreeWrapper]:
    """
    Run the analysis in consensus mode.

    Required more than 1 files. Return the tree and wrapper for rerun purpose.
    """
    wrapper = MFA2TreeWrapper(inputs, seqtype=seqtype)
    logger.info(f"Use {_config.avail_tree_methods[method]} to generate trees for treeness filtering.")
    wrapper.build(method=method, threads=threads)
    wrapper.compute_toverr(threads=threads)
    _output_treeness(wrapper, output, n=top_n_toverr)

    if 0 < top_n_toverr < len(wrapper):
        corrected_top_n = len(wrapper[:top_n_toverr])
        logger.info(f"Use the MSAs which the treeness are among the top {corrected_top_n} for final tree building.")
        wrapper = wrapper[:corrected_top_n]
    else:
        logger.info("Use all MSAs for final tree building.")

    logger.info("Generate trees on selected MSAs and conclude a majority consensus tree.")
    return wrapper.get_consensus_tree(), wrapper


def determine_samples_and_seqtype(inputs_dir: Path | None) -> tuple[tuple[str] | None, str | None]:
    """
    Determine the sample names and seqtype from libphyling checkpoint.

    Autodetect the names/seqtype from sequences if checkpoint is not available.
    """
    if inputs_dir:
        libphylingPrecheck.setup(folder=inputs_dir)
        try:
            _, samplelist, _ = libphylingPrecheck.load_checkpoint()
            seqtype = samplelist.seqtype
            samples = tuple(sample.name for sample in samplelist)
            logger.info(f"Inputs are {seqtype} sequences.")
            return samples, seqtype
        except RuntimeError:
            msg = "Align module checkpoint file corrupted. "
        except FileNotFoundError:
            msg = "Align module checkpoint file not found. "
        except Exception:
            msg = "Unknown issue to read samples and seqtype info from checkpoint file."
    else:
        msg = "Have ambiguous input directory. "
    logger.warning(msg + "Determine samples and seqtype from the given inputs.")
    return None, None


def _output_treeness(wrapper: MFA2TreeWrapper, output: Path, *, n: int = 0) -> None:
    """Output the treeness files."""
    if len(wrapper) == 1:
        mfa2tree_list = [wrapper[0]]
    else:
        mfa2tree_list: list[MFA2Tree] = wrapper[:n]
    treeness_tsv = output / expected_files["treeness"]
    with open(treeness_tsv, "w") as f:
        for mfa2tree in mfa2tree_list:
            print(mfa2tree.name, mfa2tree.toverr, sep="\t", file=f)
    logger.info(f"Selected markers and their treeness/RCV scores are output to {treeness_tsv}")

    # Symlink to seletced MSAs
    msas_dir = output / expected_files["msas_dir"]
    msas_dir.mkdir(exist_ok=True)
    inputs = [mfa2tree.path for mfa2tree in mfa2tree_list]
    for file in inputs:
        (msas_dir / file.name).symlink_to(file.absolute())
