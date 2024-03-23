"""Phylogenetic tree construction methods."""

from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import hashlib
import logging
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
import time
from io import StringIO
from itertools import product
from multiprocessing.dummy import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import NCBICodonTable, NCBICodonTableDNA
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phykit.services.alignment.base import Alignment as Phykit_alignment
from phykit.services.tree.base import Tree as Phykit_tree

import phyling.config
import phyling.error


class MFA2Tree:
    """Convert peptide/DNA multiple sequence alignment to Biopython Tree object and calculate the treeness/RCV of it."""

    def __init__(self, mfa: Path, seqtype: str, partition: Path | None = None):
        """Initialize the MFA2Tree object."""
        self._mfa = mfa
        if partition:
            self._partition = partition
        supported_seqtype = (phyling.config.seqtype_pep, phyling.config.seqtype_dna)
        if seqtype not in supported_seqtype:
            raise KeyError(f'argument "seqtype" not falls in {supported_seqtype}')
        self._seqtype = seqtype

    @property
    def mfa(self) -> Path:
        """Return the mfa path."""
        return self._mfa

    @property
    def tree(self) -> Phylo.BaseTree.Tree:
        """Return the Biopython tree object."""
        if not hasattr(self, "_tree"):
            raise AttributeError("No self._tree found. Please make sure the build function was run successfully.")
        return self._tree

    def build(self, method: str, output: Path | None = None, threads: int = 1, verbose=False) -> None:
        """Run tree building step and return the Biopython Tree object."""
        if method not in phyling.config.avail_tree_methods.keys():
            raise KeyError(f'Argument "method" not falls in {phyling.config.avail_tree_methods.keys()}')
        if method in ("upgma", "nj"):
            tree = self._with_phylo_module(method=method)
        elif method == "ft":
            tree = self._with_FastTree(verbose=verbose)
        elif method == "raxml":
            tree = self._with_RAxML(output=output, threads=threads, verbose=verbose)
        elif method == "iqtree":
            tree = self._with_IQTree(output=output, threads=threads, verbose=verbose)
        else:
            raise NotImplementedError(f"{phyling.config.avail_tree_methods[method]} is not implemented yet.")
        logging.debug(f"Tree building on {self._mfa.name} is done.")
        self._tree = tree

    def treeness_over_rcv(self) -> float:
        """Calculate the treeness/RCV of the tree by PhyKIT implementation."""
        # calculate treeness
        treeness = Phykit_tree().calculate_treeness(tree=self._tree)

        # calculate rcv
        aln = Phykit_alignment(alignment_file_path=self._mfa)
        relative_composition_variability = aln.calculate_rcv()

        # calculate treeness/rcv
        self._toverr = round(treeness / relative_composition_variability, 4)
        return self._toverr

    def _with_phylo_module(self, method) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using a simple distance method."""
        msa = AlignIO.read(self._mfa, format="fasta")
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, method)
        return constructor.build_tree(msa)

    def _with_VeryFastTree(self, threads: int) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using VeryFastTree."""
        if self._seqtype == phyling.config.seqtype_dna:
            cmd = ["VeryFastTree", "-nt", "-gamma", "-threads", str(threads), self._mfa]
        else:
            cmd = ["VeryFastTree", "-lg", "-gamma", "-threads", str(threads), self._mfa]
        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        tree, err = p.communicate()
        if not tree:
            print(err)
            sys.exit(1)
        return Phylo.read(StringIO(tree), "newick")

    def _with_FastTree(self, verbose: bool) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using FastTree."""
        fasttree_bins = ("FastTree", "fasttree")

        for fasttree in fasttree_bins:
            fasttree = shutil.which(fasttree)
            if fasttree and verbose:
                logging.debug(f"Found FastTree at {fasttree}")
                break

        if not fasttree:
            raise phyling.error.BinaryNotFoundError(
                'FastTree not found. Please install it through "conda install -c bioconda fasttree"'
            )

        cmd = [fasttree, "-gamma", str(self._mfa)]
        if self._seqtype == phyling.config.seqtype_dna:
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
            print(err)
            sys.exit(1)
        return Phylo.read(StringIO(tree), "newick")

    def _with_RAxML(self, output: Path, threads: int, verbose: bool) -> Phylo.BaseTree.Tree:
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
            raise phyling.error.BinaryNotFoundError(
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
            str(self._mfa),
            "-p",
            "12345",
            "-w",
            str((output_raxml).absolute()),
        ]
        if self._seqtype == phyling.config.seqtype_dna:
            seqtype = "cds"
            cmd.extend(["-m", "GTRCAT", "-n", seqtype])
        else:
            seqtype = "pep"
            cmd.extend(["-m", "PROTGAMMAAUTO", "-n", seqtype])
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
            tree = Phylo.read(output_raxml / f"RAxML_bestTree.{seqtype}", "newick")
        except FileNotFoundError:
            print(str(self._mfa))
            print(stdout)
            sys.exit(1)
        if output is None:
            shutil.rmtree(tempdir)

        return tree

    def _with_IQTree(self, output: Path, threads: int, verbose: bool) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using IQTree."""
        iqtree = shutil.which("iqtree")
        if iqtree and verbose:
            logging.debug(f"Found IQTree at {iqtree}")
        if not iqtree:
            raise phyling.error.BinaryNotFoundError(
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
            str(self._mfa),
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
            print(str(self._mfa))
            print(stderr)
            sys.exit(1)

        try:
            tree = Phylo.read(output_iqtree / "iqtree.treefile", "newick")
        except FileNotFoundError:
            print(str(self._mfa))
            print(stdout)
            sys.exit(1)

        if output is None:
            shutil.rmtree(tempdir)

        return tree


class SingleTreeGenerator:
    """A tree generator that for single alignment input."""

    def __init__(self, mfa: Path, seqtype: str, partition: Path | None = None):
        """Initialize the tree generator object."""
        self._seqtype = seqtype
        self._mfa2tree_obj_list = [MFA2Tree(mfa, seqtype=seqtype, partition=partition)]

    def build(self, method: str, output: Path | None = None, threads: int = 1, verbose=False) -> None:
        """Run the tree building step for each MFA2Tree object."""
        func_start = time.monotonic()
        logging.info(f"Start tree building with {phyling.config.avail_tree_methods[method]}...")

        if len(self._mfa2tree_obj_list) == 1 or threads == 1:
            logging.debug(f"Run in sequential mode with maximum {threads} threads.")
            self._trees = []
            for mfa2tree_obj in self._mfa2tree_obj_list:
                self._trees.append(
                    self._build_helper(mfa2tree_obj, method=method, output=output, threads=threads, verbose=verbose)
                )
        else:
            logging.debug(f"Run in multiprocesses mode. {threads} jobs are run concurrently.")
            with Pool(threads) as pool:
                self._trees = pool.starmap(
                    self._build_helper,
                    ((mfa2tree_obj, method, output, threads, verbose) for mfa2tree_obj in self._mfa2tree_obj_list),
                )
        logging.debug(
            f"Tree building with {phyling.config.avail_tree_methods[method]} was finished in "
            f"{phyling.config.runtime(func_start)}."
        )

    def get_tree(self) -> Phylo.BaseTree.Tree:
        """Get the tree."""
        return self._mfa2tree_obj_list[0].tree

    def get_mfa(self) -> Path:
        """Get the mfa path."""
        return self._mfa2tree_obj_list[0].mfa

    def _build_helper(
        self, mfa2tree_obj: MFA2Tree, method: str, output: Path | None, threads: int, verbose: bool
    ) -> Phylo.BaseTree.Tree:
        """Helper function to run the MFA2tree build method."""
        return mfa2tree_obj.build(method=method, output=output, threads=threads, verbose=verbose)


class TreesGenerator(SingleTreeGenerator):
    """A trees generator that support treeness/RCV estimation, consensus tree and MSAs concatenation feature."""

    def __init__(self, *mfas: Path, seqtype: str):
        """Initialize the tree generator object."""
        self._seqtype = seqtype
        self._mfa2tree_obj_list = [MFA2Tree(mfa, seqtype=seqtype) for mfa in mfas]

    def treeness_over_rcv(self, threads: int = 1) -> None:
        """Calculate the treeness/RCV for each MFA2Tree object and sort the MFA2Tree object by the values."""
        func_start = time.monotonic()
        logging.info("Calculating treeness/RCV...")
        if len(self._mfa2tree_obj_list) == 1 or threads == 1:
            logging.debug("Run in sequential mode.")
            self._toverr = []
            for mfa2tree_obj in self._mfa2tree_obj_list:
                self._toverr.append(self._toverr_helper(mfa2tree_obj))
        else:
            logging.debug(f"Run in multiprocesses mode. {threads} jobs are run concurrently.")
            with Pool(threads) as pool:
                self._toverr = pool.map(self._toverr_helper, self._mfa2tree_obj_list)
        self._sort_by_toverr()
        logging.debug(f"Treeness/RCV estimation was finished in {phyling.config.runtime(func_start)}.")

    def get_tree(self, n: int = 0) -> Phylo.BaseTree.Tree | list[Phylo.BaseTree.Tree]:
        """Get the top n trees based on treeness/RCV."""
        n = self._get_helper(n)
        if n == 1:
            return self._mfa2tree_obj_list[0].tree
        return [mfa2tree_obj.tree for mfa2tree_obj in self._mfa2tree_obj_list[:n]]

    def get_mfa(self, n: int = 0) -> Path | list[Path]:
        """Get the mfa paths of the top n trees based on treeness/RCV."""
        n = self._get_helper(n)
        if n == 1:
            return self._mfa2tree_obj_list[0].mfa
        return [mfa2tree_obj.mfa for mfa2tree_obj in self._mfa2tree_obj_list[:n]]

    def get_toverr(self, n: int = 0) -> float | list[float]:
        """Get the top n treeness/RCVs."""
        n = self._get_helper(n)
        if n == 1:
            return self._toverr[0]
        return self._toverr[:n]

    def concat(self, output: Path, taxonList: list, n: int = 0, partition: str = None, threads: int = 1) -> None:
        """Concatenate selected MSAs to a single mfa file and return the corresponding SingleTreeGenerator object."""
        func_start = time.monotonic()
        mfas = self.get_mfa(n=n)
        if not isinstance(mfas, list):
            raise ValueError("Only one MSA found. Don't need to run concat.")
        logging.info("Concatenate selected MSAs...")
        alignmentList = []
        for mfa in mfas:
            alignment = AlignIO.read(mfa, format="fasta")
            # Assign model for raxml/iqtree compatible partition file
            alignment.annotations["model"] = "DNA" if self._seqtype == phyling.config.seqtype_dna else "AUTO"
            alignment.annotations["seqname"] = mfa.name
            alignmentList.append(alignment)
        # Sort by MSA length, since shorter MSA might have higher chances to have missing states
        alignmentList = sorted(alignmentList, key=lambda alignment: alignment.get_alignment_length())

        # Create an empty MSA holder
        concat_alignments = MultipleSeqAlignment([])
        for sample in taxonList:
            concat_alignments.append(SeqRecord(Seq(""), id=sample, description=""))
        concat_alignments.sort()

        # Fill missing taxon on some of the MSAs
        if threads > 1:
            with Pool(threads) as pool:
                alignmentList = pool.starmap(self._fill_missing_taxon, product([taxonList], alignmentList))
        else:
            alignmentList = [self._fill_missing_taxon(taxonList, alignment) for alignment in alignmentList]
        logging.debug("Filling missing taxon done.")

        # Save partition_info
        if partition:
            alignmentList = self._concat_msa_with_missing_states(alignmentList)
            logging.debug("Missing states checking done.")
            if len(alignmentList) == 1:
                logging.error("Only 1 MSA left after combining msas with missing states. Please rerun with --partition disabled.")
                sys.exit(1)

        # Concatenate
        for alignment in alignmentList:
            concat_alignments += alignment
        if partition:
            concat_alignments.annotations["partition"] = self._get_partition_info(alignmentList, partition)

        for seq in concat_alignments:
            seq.description = ""
        logging.debug("Concatenate fasta done.")

        concat_file, partition_file = self._output_concat_file(output, concat_alignments)
        logging.debug(f"Fasta concatenation was finished in {phyling.config.runtime(func_start)}.")

        return SingleTreeGenerator(concat_file, partition=partition_file, seqtype=self._seqtype)

    def get_consensus_tree(self, n: int = 0) -> Phylo.BaseTree.Tree:
        """Run astral to get consensus tree."""
        func_start = time.monotonic()
        trees = self.get_tree(n=n)
        if not isinstance(trees, list):
            logging.warning("Only found 1 tree. Don't need to run ASTRAL.")
            return trees

        if not shutil.which("astral"):
            raise phyling.BinaryNotFoundError(
                "ASTRAL not found. "
                "Please build the C++ version from the source following the instruction on https://github.com/chaoszhang/ASTER"
            )

        logging.info("Run ASTRAL to resolve consensus among multiple trees.")
        temp = StringIO()
        Phylo.write(trees, temp, "newick")
        temp.seek(0)
        p = subprocess.Popen(
            ["astral", "/dev/stdin"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, _ = p.communicate(temp.read())
        temp.close()
        logging.debug(f"Consensus tree estimation with ASTRAL was finished in {phyling.config.runtime(func_start)}.")
        return Phylo.read(StringIO(stdout), "newick")

    def _sort_by_toverr(self) -> None:
        """Sort the trees/mfas by treeness/RCV."""
        sorted_idx = np.argsort(self._toverr)[::-1]
        self._toverr = [self._toverr[idx] for idx in sorted_idx]
        self._mfa2tree_obj_list = [self._mfa2tree_obj_list[idx] for idx in sorted_idx]

    def _toverr_helper(self, mfa2tree_obj: MFA2Tree) -> float:
        """Helper function to run the MFA2tree treeness_over_rcv method."""
        return mfa2tree_obj.treeness_over_rcv()

    def _get_helper(self, n: int) -> int:
        """Validate whether the given n is within the range."""
        if n:
            if not hasattr(self, "_toverr"):
                raise AttributeError("Need to run the function treeness_over_rcv first when specifying argument n.")
        if not 0 < n <= len(self._mfa2tree_obj_list):
            if n == 0:
                logging.debug("By default use all trees.")
            else:
                logging.debug("Argument --top_n_toverr/-n is more than available mfas. Use all mfas.")
            n = len(self._mfa2tree_obj_list)
        return n

    def _fill_missing_taxon(self, taxonList: list, alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
        """Include empty gap-only strings in alignment for taxa lacking an ortholog."""
        missing = set(taxonList) - {seq.id for seq in alignment}
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

    def _concat_msa_with_missing_states(self, alignmentList: list[MultipleSeqAlignment]) -> list[MultipleSeqAlignment]:
        """Concatenate the msas with missing states to prevent warning from RAxML."""
        model = alignmentList[0].annotations["model"]

        idx = 0
        count = 0
        alignmentList_length_orig = len(alignmentList)
        makeup_partition = MultipleSeqAlignment([], annotations={"model": model})
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
        chars = set()
        for rec in alignment:
            chars.update(set(rec.seq))

        if "-" in chars:
            chars.remove("-")

        All_chars = (
            set(NCBICodonTable.protein_alphabet)
            if alignment.annotations["model"] == "AUTO"
            else set(NCBICodonTableDNA.nucleotide_alphabet)
        )
        missing_chars = All_chars - chars

        return missing_chars if missing_chars else None

    def _get_partition_info(self, alignmentList: list[MultipleSeqAlignment], partition: str):
        end = 0
        partition_info = []
        model = alignmentList[0].annotations["model"]

        if partition in ("codon", "seq+codon") and model == "AUTO":
            logging.warning('Peptide sequence detected. Force to use "seq" mode for partitioning.')
            partition = "seq"
        elif partition:
            logging.info(f"Partitioning with {partition} mode.")

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
        concat_file = output / f"{phyling.config.concat_file_basename}.{phyling.config.aln_ext}"
        with open(concat_file, "w") as f:
            SeqIO.write(concat_alignments, f, format="fasta")
            logging.info(f"Concatenated fasta is output to {concat_file}")
        if "partition" in concat_alignments.annotations:
            partition_file = output / f"{phyling.config.concat_file_basename}.{phyling.config.partition_ext}"
            with open(partition_file, "w") as f:
                [print(line, file=f) for line in concat_alignments.annotations["partition"]]
            logging.info(f"Partition file is output to {partition_file}")
        else:
            partition_file = None
            logging.debug("The concat_alignments object doesn't have partition info. No partition file output.")
        return concat_file, partition_file


def mfa_to_finaltree(
    inputs: list,
    output: Path,
    method: str,
    seqtype: str,
    samples: list,
    concat: bool = False,
    top_n_toverr: int = 50,
    partition: str | None = None,
    threads: int = 1,
    rerun: int = 0,
) -> tuple[Phylo.BaseTree.Tree, TreesGenerator]:
    """Pipeline that generate the final single tree (by either concatenate or consensus methods) from MFAs."""
    logging.info("Generate phylogenetic tree on MSA fasta.")

    if rerun == 0:
        tree_generator_obj = TreesGenerator(*inputs, seqtype=seqtype)
        if len(inputs) == 1:
            tree_generator_obj.build(method=method, output=output, threads=threads, verbose=True)
            top_n_toverr = 1
        if concat:
            logging.info("Use FastTree to generate trees for treeness filtering.")
            tree_generator_obj.build(method="ft", threads=threads)
        else:
            logging.info(f"Use {phyling.config.avail_tree_methods[method]} to generate trees for treeness filtering.")
            tree_generator_obj.build(method=method, threads=threads)

        tree_generator_obj.treeness_over_rcv(threads=threads)
    else:
        _, tree_generator_obj = _load_phylotree_ckp(output)

    if rerun <= 1:
        # Output treeness.tsv
        toverrs = tree_generator_obj.get_toverr(n=top_n_toverr)
        treeness_tsv = output / phyling.config.treeness_file
        with open(treeness_tsv, "w") as f:
            for input, toverr in zip(inputs, toverrs):
                print(input, toverr, sep="\t", file=f)
        logging.info(f"Selected markers and their treeness/RCV scores are output to {treeness_tsv}")

        # Symlink to seletced MSAs
        inputs = tree_generator_obj.get_mfa(n=top_n_toverr)
        selected_msas_dir = output / phyling.config.selected_msas_dir
        selected_msas_dir.mkdir(exist_ok=True)
        for file in inputs:
            (selected_msas_dir / file.name).symlink_to(file.absolute())

    if top_n_toverr == 0:
        logging.info("Use all MSAs for final tree building.")
    elif top_n_toverr == 1:
        logging.info("Only one MSA is selected. Report the tree directly.")
        return tree_generator_obj.get_tree(n=top_n_toverr)
    else:
        corrected_top_n = len(tree_generator_obj.get_tree(n=top_n_toverr))
        logging.info(f"Use the MSAs which the treeness are among the top {corrected_top_n} for final tree building.")

    if concat:
        if rerun <= 2:
            concat_tree_obj = tree_generator_obj.concat(
                output=output, taxonList=samples, n=top_n_toverr, partition=partition, threads=threads
            )
        else:
            partition_file = (
                output / f"{phyling.config.concat_file_basename}.{phyling.config.partition_ext}" if partition else None
            )
            concat_tree_obj = SingleTreeGenerator(
                output / f"{phyling.config.concat_file_basename}.{phyling.config.aln_ext}",
                partition=partition_file,
                seqtype=seqtype,
            )
    else:
        if partition:
            logging.warning("Concatenate option disabled. No partition file will be created.")

    if concat:
        logging.info("Use the concatednated fasta to generate final tree.")
        concat_tree_obj.build(method=method, output=output, threads=threads, verbose=True)
        final_tree = concat_tree_obj.get_tree()
    else:
        logging.info("Generate trees on selected MSAs and conclude a majority consensus tree.")
        final_tree = tree_generator_obj.get_consensus_tree(n=top_n_toverr)

    return final_tree, tree_generator_obj


def _check_inputs(inputs: list[Path] | None, input_dir: Path | None) -> tuple[list[Path], Path | None]:
    if input_dir:
        inputs = [file for file in input_dir.glob(f"*.{phyling.config.aln_ext}")]
    else:
        inputs = [file for file in inputs if str(file).endswith(phyling.config.aln_ext)]
        input_dir = {file.parent for file in inputs}

        if len(input_dir) > 1:
            logging.warning("Won't load align module checkpoint file since the inputs do not came from the same folder.")
            logging.warning(input_dir)
            input_dir = None
        else:
            input_dir = input_dir.pop()
    return inputs, input_dir


def _multiple_files_md5sum(*files):
    """Calculate md5sum of multiple files."""
    hasher = hashlib.md5()
    for file in files:
        with open(file, "rb") as f:
            hasher.update(f.read())

    return hasher.hexdigest()


def _determine_samples_and_seqtype(inputs: list[Path], inputs_dir: Path | None) -> tuple[list[str], str]:
    """Load libphyling checkpoint file."""
    if inputs_dir:
        libphyling_ckp = inputs_dir / phyling.config.libphyling_checkpoint
        if Path(libphyling_ckp).exists():
            with open(libphyling_ckp, "rb") as f:
                try:
                    samples, _ = pickle.load(f)
                    seqtype = tuple(samples.values())[0].seqtype
                    samples = list(samples.keys())
                    return samples, seqtype
                except ValueError:
                    msg = "Align module checkpoint file corrupted. "
        else:
            msg = "Align module checkpoint file not found. "
    else:
        msg = "Have ambiguous inputs directory. "
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
        seqtype = phyling.config.seqtype_dna
    elif chars.issubset(set(NCBICodonTable.protein_alphabet)):
        seqtype = phyling.config.seqtype_pep
    else:
        logging.error("Cannot determine seqtype. Aborted.")
        sys.exit(1)
    return samples, seqtype


def _load_phylotree_ckp(output: Path) -> list[dict, TreesGenerator]:
    """Load phylotree checkpoint file."""
    checkpoint = output / phyling.config.phylotree_checkpoint
    with open(checkpoint, "rb") as f:
        params, tree_generator_obj = pickle.load(f)
    return params, tree_generator_obj


def _save_phylotree_ckp(output: Path, params: dict, tree_generator_obj: TreesGenerator) -> None:
    """Save phylotree checkpoint file."""
    checkpoint = output / phyling.config.phylotree_checkpoint
    with open(checkpoint, "wb") as f:
        pickle.dump((params, tree_generator_obj), f)


def _output_precheck(output: Path, params: dict) -> int:
    """Check whether the checkpoint file is exist and determine the steps needed to be rerun."""
    rerun = 0
    if not output.exists():
        output.mkdir()
    else:
        if output.is_dir():
            if (output / phyling.config.phylotree_checkpoint).exists():
                old_mfas = [file for file in (output / phyling.config.selected_msas_dir).glob(f"*.{phyling.config.aln_ext}")]
                params[f"{phyling.config.selected_msas_dir}_md5sum"] = _multiple_files_md5sum(*old_mfas)
                old_params, _ = _load_phylotree_ckp(output)
                diff_params = {param[0] for param in set(params.items()) ^ set(old_params.items())}
                logging.info("Load previous parameters from checkpoint to determine the rerun status.")
                logging.debug(diff_params)
                if not diff_params:
                    logging.info("Files not changed and parameters are identical to the previous run. Aborted.")
                    sys.exit(0)
                if "inputs_md5sum" in diff_params:
                    logging.error(f"Inputs are changed but the output directory {output} is not empty. Aborted.")
                    sys.exit(1)
                if params["concat"]:
                    if "concat" in diff_params:
                        rerun = 0
                        logging.info("Rerun all. Remove all previous files.")
                        [shutil.rmtree(x) for x in output.iterdir() if x.is_dir()]
                        [x.unlink(missing_ok=True) for x in output.iterdir() if x.is_file()]
                    elif any(param in (f"{phyling.config.selected_msas_dir}_md5sum", "top_n_toverr") for param in diff_params):
                        rerun = 1
                        logging.info("Rerun treeness/RCV filtering.")
                        [shutil.rmtree(x) for x in output.iterdir() if x.is_dir()]
                        [
                            x.unlink(missing_ok=True)
                            for x in output.iterdir()
                            if x.is_file() and x.name != phyling.config.phylotree_checkpoint
                        ]
                    elif "partition" in diff_params:
                        rerun = 2
                        logging.info("Rerun concatenate step.")
                        [
                            shutil.rmtree(x)
                            for x in Path(output).iterdir()
                            if x.is_dir() and x.name != phyling.config.selected_msas_dir
                        ]
                        [
                            x.unlink(missing_ok=True)
                            for x in output.iterdir()
                            if x.is_file() and x.name not in (phyling.config.treeness_file, phyling.config.phylotree_checkpoint)
                        ]
                    elif "method" in diff_params:
                        logging.info("Rerun final tree building step.")
                        rerun = 3
                        [
                            shutil.rmtree(x)
                            for x in Path(output).iterdir()
                            if x.is_dir() and x.name != phyling.config.selected_msas_dir
                        ]
                        (output / phyling.config.final_tree_file).unlink(missing_ok=True)
                else:
                    rerun = 0
                    logging.info("Rerun all. Remove all previous files.")
                    [shutil.rmtree(x) for x in output.iterdir() if x.is_dir()]
                    [x.unlink(missing_ok=True) for x in output.iterdir() if x.is_file()]
            elif any(output.iterdir()):
                logging.error(f"Checkpoint file not found but the output directory {output} is not empty. Aborted.")
                sys.exit(1)
        else:
            logging.error(f"{output} is already existed but not a folder. Aborted.")
            sys.exit(1)

    return rerun


def phylotree(inputs, input_dir, output, method, top_n_toverr, concat, partition, figure, threads, **kwargs):
    """
    Construct a phylogenetic tree based on the results of multiple sequence alignment (MSA).

    By default the consensus tree method will be employed which use a 50% cutoff to represent the majority of all the
    trees. You can use the -c/--concat option to concatenate the MSA and build a single tree instead. Note that enable the
    -c/--concat option will also output a partition file that compatible to RAxML and IQTree.

    For the tree building algorithm, the UPGMA will be used by default. Users can switch to the Neighbor Joining or FastTree
    by specifying the -m/--method nj/ft.

    The align step usually report a lot of markers but many of them are uninformative or susceptible to composition bias.
    The Treeness/RCV value computed by PhyKIT is used to estimate how informative the markers are. By default the
    -n/--top_n_toverr is set to 50 to select only the top 50 markers.

    Once the tree is built, an ASCII figure representing the tree will be displayed, and a treefile in Newick format
    will be generated as output. Additionally, users can choose to obtain a matplotlib-style figure using the
    -f/--figure option.
    """
    module_start = time.monotonic()
    params = {"method": method, "top_n_toverr": top_n_toverr, "concat": concat, "partition": partition}

    inputs, input_dir = _check_inputs(inputs=inputs, input_dir=input_dir)
    params["inputs_md5sum"] = _multiple_files_md5sum(*inputs)
    logging.info(f"Found {len(inputs)} MSA fasta.")
    samples, seqtype = _determine_samples_and_seqtype(inputs, input_dir)
    logging.info(f"Inputs are {seqtype} sequences.")

    output = Path(output)
    rerun = _output_precheck(output, params)
    logging.debug(f"rerun = {rerun}")

    final_tree, tree_generator_obj = mfa_to_finaltree(
        inputs=inputs,
        output=output,
        method=method,
        seqtype=seqtype,
        samples=samples,
        concat=concat,
        top_n_toverr=top_n_toverr,
        partition=partition,
        threads=threads,
        rerun=rerun,
    )

    selected_msas = [file for file in (output / phyling.config.selected_msas_dir).glob(f"*.{phyling.config.aln_ext}")]
    params[f"{phyling.config.selected_msas_dir}_md5sum"] = _multiple_files_md5sum(*selected_msas)
    _save_phylotree_ckp(output, params, tree_generator_obj)

    Phylo.draw_ascii(final_tree)

    output_tree = output / phyling.config.final_tree_file
    logging.info(f"Output tree to {output_tree}")
    with open(output_tree, "w") as f:
        if concat:
            strategy = "concatenate"
        else:
            strategy = "consensus"
        print(
            f"# Final tree is built using {phyling.config.avail_tree_methods[method]} with {strategy} strategy",
            file=f,
        )
        if concat and partition:
            print(f"# Partition is enabled using {partition} mode", file=f)
        Phylo.write(final_tree, f, "newick")

    if figure:
        fig, ax = plt.subplots(figsize=(20, 12))
        output_fig = output / f"{method}_tree.png"
        logging.info(f"Output figure to {output_fig}")
        Phylo.draw(final_tree, axes=ax)
        fig.savefig(output_fig)

    logging.debug(f"Tree module finished in {phyling.config.runtime(module_start)}.")
