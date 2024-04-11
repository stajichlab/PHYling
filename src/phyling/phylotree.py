"""Phylogenetic tree construction methods."""

from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import hashlib
import logging
import re
import shutil
import subprocess
import tempfile
import time
from io import StringIO
from itertools import product
from multiprocessing.dummy import Pool
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import NCBICodonTable, NCBICodonTableDNA
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phykit.services.alignment.base import Alignment as Phykit_alignment
from phykit.services.tree.base import Tree as Phykit_tree

import phyling._abc as _abc
import phyling.config as config
import phyling.exception as exception
import phyling.utils as utils
from phyling.libphyling import OutputPrecheck as libphylingPrecheck

expected_files = {
    "msas_dir": "selected_MSAs",
    "treeness": "treeness.tsv",
    "concat": f"concat_alignments.{config.aln_ext}",
    "partition": "concat_alignments.partition",
    "tree_nw": "final_tree.nw",
    "tree_img": "final_tree.png",
    "raxml": "raxml",
    "iqtree": "iqtree",
}


class MFA2Tree(_abc.SeqFileWrapperABC):
    """Convert peptide/DNA multiple sequence alignment to Biopython Tree object and calculate the treeness/RCV."""

    def __init__(self, file: str | Path, seqtype: str) -> None:
        """Initialize the object and check sequence type."""
        super().__init__(file, seqtype)

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
        return self.toverr < other.toverr

    @property
    def name(self) -> str:
        """Return the basename of the file."""
        return self.path.name

    @property
    def checksum(self) -> int:
        """Compute the checksum of the stored data."""
        return utils.get_file_checksum(self.path)

    @property
    def msa(self) -> MultipleSeqAlignment:
        """Return Biopython MultipleSeqAlignment object."""
        alignment = AlignIO.read(self.path, format="fasta")
        alignment.annotations["seqtype"] = self.seqtype
        alignment.annotations["seqname"] = self.name
        return alignment

    @property
    def tree(self) -> Tree:
        """Return the Biopython tree object."""
        if not hasattr(self, "_tree"):
            raise AttributeError("build() need to be run first.")
        return self._tree

    @property
    def method(self) -> str:
        """Return the treeness/RCVs."""
        if not hasattr(self, "_method"):
            raise AttributeError("build() need to be run first.")
        return self._method

    @property
    def toverr(self) -> float:
        """Return the treeness/RCVs."""
        if not hasattr(self, "_toverr"):
            raise AttributeError("compute_toverr() need to be run first.")
        return self._toverr

    def build(self, method: str, output: Path | None = None, threads: int = 1, verbose=False) -> None:
        """Run tree building step and return the Biopython Tree object."""
        if method not in config.avail_tree_methods.keys():
            raise KeyError(f'Argument "method" not falling in {config.avail_tree_methods.keys()}')
        self._method = method
        if method in ("upgma", "nj"):
            tree = self._with_phylo_module(method=method)
        elif method == "ft":
            tree = self._with_FastTree(verbose=verbose)
        elif method == "raxml":
            tree = self._with_RAxML(output=output, threads=threads, verbose=verbose)
        elif method == "iqtree":
            tree = self._with_IQTree(output=output, threads=threads, verbose=verbose)
        else:
            raise NotImplementedError(f"{config.avail_tree_methods[method]} is not implemented yet.")
        logging.debug(f"Tree building on {self.name} is done.")
        self._tree = tree

    def compute_toverr(self) -> None:
        """Calculate the treeness/RCV of the tree by PhyKIT implementation."""
        # calculate treeness
        treeness = Phykit_tree().calculate_treeness(tree=self.tree)

        # calculate rcv
        aln = Phykit_alignment(alignment_file_path=self.path)
        relative_composition_variability = aln.calculate_rcv()

        # calculate treeness/rcv
        self._toverr = round(treeness / relative_composition_variability, 4)

    def _with_phylo_module(self, method) -> Tree:
        """Run the tree calculation using a simple distance method."""
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, method)
        return constructor.build_tree(self.msa)

    def _with_VeryFastTree(self, threads: int) -> Tree:
        """Run the tree calculation using VeryFastTree."""
        if self._seqtype == config.seqtype_cds:
            cmd = ["VeryFastTree", "-nt", "-gamma", "-threads", str(threads), self.path]
        else:
            cmd = ["VeryFastTree", "-lg", "-gamma", "-threads", str(threads), self.path]
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

    def _with_FastTree(self, verbose: bool) -> Tree:
        """Run the tree calculation using FastTree."""
        fasttree_bins = ("FastTree", "fasttree")

        for fasttree in fasttree_bins:
            fasttree = shutil.which(fasttree)
            if fasttree and verbose:
                logging.debug(f"Found FastTree at {fasttree}")
                break

        if not fasttree:
            raise exception.BinaryNotFoundError(
                'FastTree not found. Please install it through "conda install -c bioconda fasttree"'
            )

        cmd = [fasttree, "-gamma", str(self.path)]
        if self.seqtype == config.seqtype_cds:
            cmd.insert(1, "-nt")
        else:
            cmd.insert(1, "-lg")
        if verbose:
            logging.info(f'Tree building cmd: {" ".join(cmd)}')

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

    def _with_RAxML(self, output: Path, threads: int, verbose: bool) -> Tree:
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
            if raxml and verbose:
                logging.debug(f"Found RAxML at {raxml}")
                break

        if not raxml:
            raise exception.BinaryNotFoundError(
                'RAxML not found. Please install it through "conda install -c bioconda raxml"'
            )

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
        if self.seqtype == config.seqtype_cds:
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
            logging.info(f'Tree building cmd: {" ".join(cmd)}')

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

    def _with_IQTree(self, output: Path, threads: int, verbose: bool) -> Tree:
        """Run the tree calculation using IQTree."""
        iqtree = shutil.which("iqtree")
        if iqtree and verbose:
            logging.debug(f"Found IQTree at {iqtree}")
        if not iqtree:
            raise exception.BinaryNotFoundError(
                'IQTree not found. Please install it through "conda install -c bioconda iqtree"'
            )

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
            logging.info(f'Tree building cmd: {" ".join(cmd)}')

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


class ConcatMFA2Tree(MFA2Tree):
    """
    Convert concatenated peptide/DNA multiple sequence alignment to Biopython Tree object and calculate the treeness/RCV.

    Optionally take a partition file to build tree with partition mode with RAxML and IQTree.
    """

    def __init__(self, file: str | Path, seqtype: str, partition: Path | None = None) -> None:
        """Initialize the object and check sequence type."""
        super().__init__(file, seqtype)
        if partition:
            self._partition = partition


class MFA2TreeWrapper(_abc.DataListABC[MFA2Tree]):
    """A Wrapper for MFA2Tree object(s) which also provides concatenated data preparation and consensus tree calculation."""

    def __init__(self, files: list[str | Path | MFA2Tree], seqtype: str) -> None:
        """Initialize the wrapper object."""
        super().__init__(files)
        for file in files:
            if isinstance(file, MFA2Tree):
                self.data.append(file)
            elif isinstance(file, str | Path):
                self.data.append(MFA2Tree(file, seqtype=seqtype))
            else:
                raise TypeError(f"{self.__class__.__qualname__} only accepts list of str/Path/MFA2Tree")
        self._seqtype = seqtype
        self._issorted = False

    def __repr__(self) -> str:
        """Return the string representation."""
        info = f"seqtype={self.seqtype}, sorted={self.is_sorted}"
        if self.is_sorted:
            info += f", method={self.method}"
        return super().__repr__() % info

    def __getitem__(self, idx: int | slice) -> MFA2Tree | MFA2TreeWrapper[MFA2Tree]:
        """Get the data by index or slice. Note that the data need to be sorted before getting with slice."""
        data = self.data[idx]
        if isinstance(data, MFA2Tree):
            return data
        elif isinstance(data, list):
            self._unsorted_check()
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

    def build(self, method: str, output: Path | None = None, threads: int = 1, verbose=False) -> None:
        """Run the tree building step for each MFA2Tree object."""
        func_start = time.monotonic()
        logging.info(f"Start tree building with {config.avail_tree_methods[method]}...")

        if len(self) == 1 or threads == 1:
            logging.debug(f"Sequential mode with maximum {threads} threads.")
            for mfa2tree in self:
                self._build_helper(mfa2tree, method=method, output=output, threads=threads, verbose=verbose)
        else:
            logging.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
            with Pool(threads) as pool:
                pool.starmap(
                    self._build_helper,
                    ((mfa2tree, method, output, threads, verbose) for mfa2tree in self),
                )
        logging.debug(
            f"Tree building with {config.avail_tree_methods[method]} was finished in " f"{utils.runtime(func_start)}."
        )

    def compute_toverr(self, threads: int = 1) -> None:
        """Calculate the treeness/RCV for each MFA2Tree object and sort the MFA2Tree object by the values."""
        func_start = time.monotonic()
        logging.info("Calculating treeness/RCV...")
        if len(self) == 1 or threads == 1:
            logging.debug(f"Sequential mode with maximum {threads} threads.")
            for mfa2tree in self:
                self._compute_toverr_helper(mfa2tree)
        else:
            logging.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
            with Pool(threads) as pool:
                pool.map(self._compute_toverr_helper, self)
        self.sort()
        logging.info("Done.")
        logging.debug(f"Treeness/RCV estimation was finished in {utils.runtime(func_start)}.")

    def get_consensus_tree(self) -> Tree:
        """Run astral to get the consensus tree."""
        func_start = time.monotonic()
        self._single_file_check(self.concat.__name__)
        self._unsorted_check()
        if not shutil.which("astral"):
            raise exception.BinaryNotFoundError(
                "ASTRAL not found. "
                "Please build the C++ version from the source following the instruction on https://github.com/chaoszhang/ASTER"
            )

        logging.info("Run ASTRAL to resolve consensus among multiple trees.")
        temp = StringIO()
        Phylo.write([mfa2tree.tree for mfa2tree in self], temp, "newick")
        temp.seek(0)
        p = subprocess.Popen(
            ["astral", "/dev/stdin"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, _ = p.communicate(temp.read())
        temp.close()
        logging.info("Done.")
        logging.debug(f"Consensus tree estimation with ASTRAL was finished in {utils.runtime(func_start)}.")
        return Phylo.read(StringIO(stdout), "newick")

    def concat(
        self, output: Path, samples: Iterable[str], partition: str | None = None, threads: int = 1
    ) -> ConcatMFA2Tree:
        """Concatenate selected MSAs to a single mfa file and return the corresponding ConcatMFA2Tree object."""
        func_start = time.monotonic()
        self._single_file_check(self.concat.__name__)
        self._unsorted_check()
        if self.method not in ("raxml", "iqtree") and partition:
            logging.warning("Partition is forced to be disabled since it only works when using raxml and iqtree.")
            partition = None

        alignmentList = [mfa2tree.msa for mfa2tree in self]
        # Sort by MSA length, since shorter MSA might have higher chances to have missing states
        alignmentList = sorted(alignmentList, key=lambda alignment: alignment.get_alignment_length())

        # Create an empty MSA holder
        concat_alignments = MultipleSeqAlignment([])
        for sample in samples:
            concat_alignments.append(SeqRecord(Seq(""), id=sample, description=""))
        concat_alignments.sort()

        # Fill missing taxon on some of the MSAs
        logging.info("Fill missing taxon...")
        if threads > 1:
            with Pool(threads) as pool:
                alignmentList = pool.starmap(self.fill_missing_taxon, product([samples], alignmentList))
        else:
            alignmentList = [self.fill_missing_taxon(samples, alignment) for alignment in alignmentList]
        logging.info("Done.")

        # Save partition_info
        if partition:
            logging.info("Check msa with missing states...")
            alignmentList = self._concat_msa_with_missing_states(alignmentList)
            logging.info("Done.")
            if len(alignmentList) == 1:
                raise RuntimeError(
                    "Only 1 MSA left after combining msas with missing states. Please rerun with --partition disabled."
                )

        # Concatenate
        logging.info("Concatenate selected MSAs...")
        for alignment in alignmentList:
            concat_alignments += alignment
        if partition:
            concat_alignments.annotations["partition"] = self._get_partition_info(alignmentList, partition)
        for seq in concat_alignments:
            seq.description = ""
        logging.info("Done.")

        concat_file, partition_file = self._output_concat_file(output, concat_alignments)
        logging.debug(f"Fasta concatenation was finished in {utils.runtime(func_start)}.")

        return ConcatMFA2Tree(concat_file, partition=partition_file, seqtype=self.seqtype)

    @staticmethod
    def fill_missing_taxon(samples: Iterable[str], alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
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

    def _build_helper(self, mfa2tree: MFA2Tree, method: str, output: Path | None, threads: int, verbose: bool) -> Tree:
        """Helper function to run the MFA2tree build method."""
        mfa2tree.build(method=method, output=output, threads=threads, verbose=verbose)

    def _compute_toverr_helper(self, mfa2tree: MFA2Tree) -> None:
        """Helper function to run the MFA2Tree compute_toverr method."""
        mfa2tree.compute_toverr()

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
            logging.warning(msg)
            logging.warning(
                "This indicates samples might contain many short sequences. Recommended to run with --partition disabled."
            )
        else:
            logging.info(msg)
        logging.info("MSAs with missing states would be concatenated into a single partition.")

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
            if self.seqtype == config.seqtype_pep
            else set(NCBICodonTableDNA.nucleotide_alphabet)
        )
        missing_chars = All_chars - chars

        return missing_chars if missing_chars else None

    def _get_partition_info(self, alignmentList: list[MultipleSeqAlignment], partition: str):
        """Get raxml-supported partition info."""
        end = 0
        partition_info = []

        if partition in ("codon", "seq+codon") and self.seqtype == config.seqtype_cds:
            logging.warning('Peptide sequence detected. Force to use "seq" mode for partitioning.')
            partition = "seq"
        if partition:
            logging.info(f"Partitioning with {partition} mode.")

        # Assign model for raxml/iqtree compatible partition file
        model = "DNA" if self.seqtype == config.seqtype_cds else "AUTO"
        for alignment in alignmentList:
            start = end + 1
            end += alignment.get_alignment_length()
            if partition == "seq":
                partition_info.append(f'{model}, {alignment.annotations["seqname"]} = {start}-{end}')
            elif partition == "seq+codon":
                for pos in range(3):
                    partition_info.append(
                        f'{model}, {alignment.annotations["seqname"]}_codon{pos + 1} = {start + pos}-{end}\\3'
                    )
        if partition == "codon":
            for pos in range(3):
                partition_info.append(f"{model}, codon{pos + 1} = {pos + 1}-{end}\\3")
        return partition_info

    def _output_concat_file(self, output: Path, concat_alignments: MultipleSeqAlignment) -> tuple[Path, Path]:
        """Output the concatenated fasta and partition file."""
        concat_file = output / expected_files["concat"]
        with open(concat_file, "w") as f:
            SeqIO.write(concat_alignments, f, format="fasta")
            logging.info(f"Concatenated fasta is output to {concat_file}")
        if "partition" in concat_alignments.annotations:
            partition_file = output / expected_files["partition"]
            with open(partition_file, "w") as f:
                [print(line, file=f) for line in concat_alignments.annotations["partition"]]
            logging.info(f"Partition file is output to {partition_file}")
        else:
            partition_file = None
            logging.debug("The concat_alignments object doesn't have partition info. No partition file output.")
        return concat_file, partition_file

    def _unsorted_check(self):
        """Check whether data is sorted and raise an error if not."""
        if not self.is_sorted:
            raise RuntimeError("Data is not sorted. Please run compute_toverr() to sort it by toverr first.")

    def _single_file_check(self, func_name: str):
        """Check whether data has more than one MFA2Tree object and raise an error if not."""
        if len(self) < 2:
            raise RuntimeError(f"{func_name} cannot be run since only single mfa found.")


class OutputPrecheck(_abc.OutputPrecheck[MFA2TreeWrapper]):
    """A class that provides features for input/output precheck, checkpoint loading/saving and final tree output."""

    def __init__(self, output: str | Path, params: dict | None = None) -> None:
        """Initialize the object with Nonetype slots for parameters and data if not receiving params and samplelist."""
        param_keys = (
            "method",
            "top_n_toverr",
            "concat",
            "partition",
            "inputs_checksum",
            f'{expected_files["msas_dir"]}_checksum',
        )
        data = None
        super().__init__(output, param_keys, data)
        self._ckp = self._output / config.phylotree_checkpoint
        if params:
            self.params = params

    @_abc.OutputPrecheck.data.setter
    def data(self, cur_data: MFA2TreeWrapper) -> None:
        """Update the current data. Only accepts MFA2TreeWrapper."""
        if not isinstance(cur_data, MFA2TreeWrapper):
            raise ValueError(f"Invalid data. Only accept {MFA2TreeWrapper} object.")
        self._data = cur_data

    def load_checkpoint(self) -> tuple[dict, MFA2TreeWrapper]:
        """
        Load the checkpoint and retrieve the required params/data to determine the rerun status.

        This should be run before precheck.
        """
        return super().load_checkpoint()

    def save_checkpoint(self) -> None:
        """Save the parameters and data as a checkpoint for rerun."""
        self.params[f'{expected_files["msas_dir"]}_checksum'] = utils.get_multifiles_checksum(
            (self.path / expected_files["msas_dir"]).glob(f"*.{config.aln_ext}")
        )
        super().save_checkpoint()

    def output_results(self, tree: Tree, figure: bool = False):
        """Output the tree in newick format and figure."""
        Phylo.draw_ascii(tree)
        output_tree = self.path / expected_files["tree_nw"]
        logging.info(f"Output tree to {output_tree}")
        with open(output_tree, "w") as f:
            strategy = "concatenate" if self.params["concat"] else "consensus"
            print(
                f'# Final tree is built using {config.avail_tree_methods[self.params["method"]]} with {strategy} strategy',
                file=f,
            )
            if self.params["concat"] and self.params["partition"]:
                print(f'# Partition is enabled using {self.params["partition"]} mode', file=f)
            Phylo.write(tree, f, "newick")

        if figure:
            fig, ax = plt.subplots(figsize=(20, 12))
            output_fig = self.path / expected_files["tree_img"]
            logging.info(f"Output figure to {output_fig}")
            Phylo.draw(tree, axes=ax)
            fig.savefig(output_fig)

    def _determine_rerun(self, prev_params: dict, prev_data: MFA2TreeWrapper) -> int:
        """Define the actions that need to do when found a checkpoint file in the given output folder."""
        super()._determine_rerun(prev_params)

        prev_mfas = [file for file in (self.path / expected_files["msas_dir"]).glob(f"*.{config.aln_ext}")]
        try:
            prev_mfas_checksum = utils.get_multifiles_checksum(prev_mfas)
        except FileNotFoundError:
            prev_mfas_checksum = -1
        cur_params = {**self.params, **{f'{expected_files["msas_dir"]}_checksum': prev_mfas_checksum}}
        diff_params = {param[0] for param in set(cur_params.items()) ^ set(prev_params.items())}
        if not diff_params:
            raise SystemExit("Files not changed and parameters are identical to the previous run. Aborted.")
        logging.debug(f"{self._determine_rerun.__name__}.diff_params = {diff_params}")
        if "inputs_checksum" in diff_params:
            raise SystemExit(f"Inputs are changed but the output directory {self.path} is not empty. Aborted.")
        if not cur_params["concat"]:
            print("Rerun all. Remove all previous files.")
            rerun_status, rm_dirs, rm_files = self._rerun_status_0()
        else:
            if "concat" in diff_params:
                logging.info("Rerun all. Remove all previous files.")
                rerun_status, rm_dirs, rm_files = self._rerun_status_0()
            elif diff_params.intersection((f'{expected_files["msas_dir"]}_checksum', "top_n_toverr")):
                logging.info("Rerun treeness/RCV filtering.")
                rerun_status, rm_dirs, rm_files = self._rerun_status_1()
            elif "partition" in diff_params:
                logging.info("Rerun concatenate step.")
                rerun_status, rm_dirs, rm_files = self._rerun_status_2()
            elif "method" in diff_params:
                logging.info("Rerun final tree building step.")
                rerun_status, rm_dirs, rm_files = self._rerun_status_3()

        if rerun_status:
            self.data = prev_data
        self._remove_files(files=rm_files, dirs=rm_dirs)
        return rerun_status

    def _rerun_status_3(self, rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
        """Return rerun status as 3 and the corresponding folders and files need to be removed."""
        for x in self.path.iterdir():
            if x.is_dir():
                if x.name != expected_files["msas_dir"]:
                    rm_dirs.add(x)
            else:
                if x.name in (expected_files["tree_nw"], expected_files["tree_img"]):
                    rm_files.add(self.path / x)
        return 3, rm_dirs, rm_files

    def _rerun_status_2(self, rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
        """Return rerun status as 2 and the corresponding folders and files need to be removed."""
        _, rm_dirs, rm_files = self._rerun_status_3(rm_dirs, rm_files)
        rm_files.update(
            x
            for x in self.path.iterdir()
            if x.is_file() and x.name not in (expected_files["treeness"], config.phylotree_checkpoint)
        )

        return 2, rm_dirs, rm_files

    def _rerun_status_1(self, rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
        """Return rerun status as 1 and the corresponding folders and files need to be removed."""
        _, rm_dirs, rm_files = self._rerun_status_2(rm_dirs, rm_files)
        for x in self.path.iterdir():
            if x.is_dir():
                rm_dirs.add(x)
            else:
                if x.name != config.phylotree_checkpoint:
                    rm_files.add(x)
        return 1, rm_dirs, rm_files

    def _rerun_status_0(self, rm_dirs: set = set(), rm_files: set = set()) -> tuple[int, set, set]:
        """Return rerun status as 0 and the corresponding folders and files need to be removed."""
        _, rm_dirs, rm_files = self._rerun_status_1(rm_dirs, rm_files)
        for x in self.path.iterdir():
            if x.is_dir():
                rm_dirs.add(x)
            else:
                rm_files.add(x)

        return 0, rm_dirs, rm_files


def single_mfa(
    input: Path,
    output: Path,
    method: str,
    seqtype: str,
    threads: int = 1,
    rerun: int = 0,
    wrapper: MFA2TreeWrapper | None = None,
) -> tuple[Tree, MFA2TreeWrapper]:
    """
    Run the analysis in single mfa mode.

    Return the tree and wrapper for rerun purpose.
    """
    if bool(rerun) ^ bool(wrapper):
        raise RuntimeError("Rerun must to be coupled with wrapper.")

    if rerun == 0:
        wrapper = MFA2TreeWrapper(input, seqtype=seqtype)
        wrapper.build(method=method, output=output, threads=threads, verbose=True)

    if rerun <= 1:
        _output_treeness(wrapper, output, threads=1)

    logging.info("Only one MSA is selected. Report the tree directly.")
    return wrapper[0].tree, wrapper


def concat(
    inputs: list,
    output: Path,
    method: str,
    seqtype: str,
    samples: list,
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
        logging.info("Use FastTree to generate trees for treeness filtering.")
        wrapper.build(method="ft", threads=threads)
        wrapper.compute_toverr(threads=threads)
    else:
        wrapper = wrapper

    if rerun <= 1:
        _output_treeness(wrapper, output, n=top_n_toverr, threads=1)

    if 0 < top_n_toverr < len(wrapper):
        corrected_top_n = len(wrapper[:top_n_toverr])
        logging.info(f"Use the MSAs which the treeness are among the top {corrected_top_n} for final tree building.")
        wrapper = wrapper[:corrected_top_n]
    else:
        logging.info("Use all MSAs for final tree building.")

    if rerun <= 2:
        concat_mfa2tree = wrapper.concat(output=output, samples=samples, partition=partition, threads=threads)
    else:
        partition_file = output / expected_files["partition"] if partition else None
        concat_mfa2tree = ConcatMFA2Tree(output / expected_files["concat"], partition=partition_file, seqtype=seqtype)

    logging.info("Use the concatednated fasta to generate final tree.")
    concat_mfa2tree.build(method=method, output=output, threads=threads, verbose=True)
    return concat_mfa2tree.tree, wrapper


def consensus(
    inputs: list,
    output: Path,
    method: str,
    seqtype: str,
    top_n_toverr: int = 50,
    threads: int = 1,
    rerun: int = 0,
    wrapper: MFA2TreeWrapper | None = None,
) -> tuple[Tree, MFA2TreeWrapper]:
    """
    Run the analysis in consensus mode.

    Required more than 1 files. Return the tree and wrapper for rerun purpose.
    """
    if bool(rerun) ^ bool(wrapper):
        raise RuntimeError("Rerun must to be coupled with wrapper.")

    if rerun == 0:
        wrapper = MFA2TreeWrapper(inputs, seqtype=seqtype)
        logging.info(f"Use {config.avail_tree_methods[method]} to generate trees for treeness filtering.")
        wrapper.build(method=method, threads=threads)
        wrapper.compute_toverr(threads=threads)
    else:
        wrapper = wrapper

    if rerun <= 1:
        _output_treeness(wrapper, output, n=top_n_toverr, threads=1)

    if 0 < top_n_toverr < len(wrapper):
        corrected_top_n = len(wrapper[:top_n_toverr])
        logging.info(f"Use the MSAs which the treeness are among the top {corrected_top_n} for final tree building.")
        wrapper = wrapper[:corrected_top_n]
    else:
        logging.info("Use all MSAs for final tree building.")

    logging.info("Generate trees on selected MSAs and conclude a majority consensus tree.")
    return wrapper.get_consensus_tree(), wrapper


def determine_samples_and_seqtype(inputs: list[Path], inputs_dir: Path | None) -> tuple[tuple[str], str]:
    """
    Determine the sample names and seqtype from libphyling checkpoint.

    Autodetect the names/seqtype from sequences if checkpoint is not available.
    """
    if inputs_dir:
        precheck_obj = libphylingPrecheck(inputs_dir)
        if precheck_obj.checkpoint.exists():
            try:
                _, samplelist, _ = precheck_obj.load_checkpoint()
                seqtype = samplelist.seqtype
                samples = tuple(sample.name for sample in samplelist)
                return samples, seqtype
            except ValueError:
                msg = "Align module checkpoint file corrupted. "
        else:
            msg = "Align module checkpoint file not found. "
    else:
        msg = "Have ambiguous input directory. "
    logging.warning(msg + "Determine samples and seqtype from the given inputs.")

    samples = set()
    chars = set()
    for file in inputs:
        alignment = AlignIO.read(file, format="fasta")
        for rec in alignment:
            samples.add(rec.id)
            chars.update(str(rec.seq))
    if "-" in chars:
        chars.remove("-")
    if len(chars) <= 4 and chars.issubset(set(NCBICodonTableDNA.nucleotide_alphabet)):
        seqtype = config.seqtype_cds
    elif chars.issubset(set(NCBICodonTable.protein_alphabet)):
        seqtype = config.seqtype_pep
    else:
        raise exception.SeqtypeError("Cannot determine seqtype. Aborted.")
    return tuple(samples), seqtype


def _output_treeness(wrapper: MFA2TreeWrapper, output: Path, n: int = 0, threads: int = 1) -> None:
    """Output the treeness files."""
    if len(wrapper) == 1:
        mfa2tree_list = [wrapper[0]]
    else:
        mfa2tree_list: list[MFA2Tree] = wrapper[:n]
    treeness_tsv = output / expected_files["treeness"]
    with open(treeness_tsv, "w") as f:
        for mfa2tree in mfa2tree_list:
            print(mfa2tree.name, mfa2tree.toverr, sep="\t", file=f)
    logging.info(f"Selected markers and their treeness/RCV scores are output to {treeness_tsv}")

    # Symlink to seletced MSAs
    msas_dir = output / expected_files["msas_dir"]
    msas_dir.mkdir(exist_ok=True)
    inputs = [mfa2tree.path for mfa2tree in mfa2tree_list]
    for file in inputs:
        (msas_dir / file.name).symlink_to(file.absolute())
