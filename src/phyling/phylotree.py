"""Phylogenetic tree construction methods."""

from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import logging
import pickle
import shutil
import subprocess
import sys
import time
from io import StringIO
from itertools import product
from multiprocessing.dummy import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phykit.services.alignment.base import Alignment as Phykit_alignment
from phykit.services.tree.base import Tree as Phykit_tree

import phyling.config
import phyling.error


class MFA2Tree:
    """Convert peptide/DNA multiple sequence alignment to Biopython Tree object and calculate the treeness/RCV of it."""

    def __init__(self, mfa: Path, seqtype: str):
        """Initialize the MFA2Tree object."""
        self._mfa = mfa
        supported_seqtype = ["peptide", "DNA"]
        if seqtype not in supported_seqtype:
            raise KeyError(f'argument "seqtype" not falls in {supported_seqtype}')
        self._seqtype = seqtype

    @property
    def mfa(self) -> Path:
        """Return the mfa path."""
        return self._mfa

    def build(self, method: str) -> None:
        """Run tree building step and return the Biopython Tree object."""
        if method not in phyling.config.avail_tree_methods.keys():
            raise KeyError(f'argument "method" not falls in {phyling.config.avail_tree_methods.keys()}')
        if method in ["upgma", "nj"]:
            tree = self._with_phylo_module(method=method)
        elif method == "ft":
            if not shutil.which("FastTree"):
                raise phyling.error.BinaryNotFoundError(
                    'FastTree not found. Please install it through "conda install -c bioconda fasttree"'
                )
            tree = self._with_FastTree(self._seqtype)
        else:
            raise NotImplementedError(f"{phyling.config.avail_tree_methods[method]} is not implemented yet.")
        logging.debug(f"Tree building on {self._mfa.name} is done")
        self._tree = tree

    @property
    def tree(self) -> Phylo.BaseTree.Tree:
        """Return the Biopython tree object."""
        if not hasattr(self, "_tree"):
            raise AttributeError("No self._tree found. Please make sure the build function was run successfully")
        return self._tree

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
        MSA = AlignIO.read(self._mfa, format="fasta")
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, method)
        return constructor.build_tree(MSA)

    def _with_VeryFastTree(self, seqtype: str, threads: int) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using VeryFastTree."""
        if seqtype == "DNA":
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

    def _with_FastTree(self, seqtype: str) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using FastTree."""
        if seqtype == "DNA":
            cmd = ["FastTree", "-nt", "-gamma", self._mfa]
        else:
            cmd = ["FastTree", "-lg", "-gamma", self._mfa]
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


class TreesGenerator:
    """A phylogentic tree generator used in phyling."""

    def __init__(self, seqtype: str, *mfas: Path):
        """Initialize the tree generator object."""
        self._mfa2tree_obj_list = [MFA2Tree(mfa, seqtype=seqtype) for mfa in mfas]

    def build(self, method: str, threads: int) -> None:
        """Run the tree building step for each MFA2Tree object."""
        logging.info("Start tree building...")
        if len(self._mfa2tree_obj_list) == 1 or threads == 1:
            logging.debug("Run in single thread mode.")
            self._trees = []
            for mfa2tree_obj in self._mfa2tree_obj_list:
                self._trees.append(self._build_helper(mfa2tree_obj, method=method))
        else:
            logging.debug(f"Run in multiprocesses mode. {threads} jobs are run concurrently")
            with Pool(threads) as pool:
                _ = pool.starmap(self._build_helper, ((mfa2tree_obj, method) for mfa2tree_obj in self._mfa2tree_obj_list))

    def treeness_over_rcv(self, threads: int) -> None:
        """Calculate the treeness/RCV for each MFA2Tree object and sort the MFA2Tree object by the values."""
        logging.info("Calculating treeness/RCV...")
        if len(self._mfa2tree_obj_list) == 1 or threads == 1:
            logging.debug("Run in single thread mode.")
            self._toverr = []
            for mfa2tree_obj in self._mfa2tree_obj_list:
                self._toverr.append(self._toverr_helper(mfa2tree_obj))
        else:
            logging.debug(f"Run in multiprocesses mode. {threads} jobs are run concurrently")
            with Pool(threads) as pool:
                self._toverr = pool.map(self._toverr_helper, self._mfa2tree_obj_list)
        self._sort_by_toverr()

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

    def _sort_by_toverr(self) -> None:
        """Sort the trees/mfas by treeness/RCV."""
        sorted_idx = np.argsort(self._toverr)[::-1]
        self._toverr = [self._toverr[idx] for idx in sorted_idx]
        self._mfa2tree_obj_list = [self._mfa2tree_obj_list[idx] for idx in sorted_idx]

    def _build_helper(self, mfa2tree_obj: MFA2Tree, method: str) -> Phylo.BaseTree.Tree:
        """Helper function to run the MFA2tree build method."""
        return mfa2tree_obj.build(method=method)

    def _toverr_helper(self, mfa2tree_obj: MFA2Tree) -> float:
        """Helper function to run the MFA2tree treeness_over_rcv method."""
        return mfa2tree_obj.treeness_over_rcv()

    def _get_helper(self, n: int) -> int:
        """Validate whether the given n is within the range."""
        if n:
            if not hasattr(self, "_toverr"):
                raise AttributeError("Need to run the function treeness_over_rcv first when specifying argument n")
        if not 0 < n <= len(self._mfa2tree_obj_list):
            if n == 0:
                logging.debug("By default use all trees")
            else:
                logging.warning("Argument --top_n_toverr/-n is larger the existing trees. Use all trees")
            n = len(self._mfa2tree_obj_list)
        return n


def mfa_to_finaltree(
    inputs: list,
    output: Path,
    method: str,
    seqtype: str,
    samples: list,
    concat: bool = False,
    top_n_toverr: int = 50,
    threads: int = 1,
) -> Phylo.BaseTree.Tree:
    """Pipeline that generate the final single tree (by either concatenate or consensus methods) from MFAs."""
    output = Path(output)
    output.mkdir(exist_ok=True)

    if top_n_toverr or len(inputs) == 1:
        func_start = time.monotonic()
        tree_generator_obj = TreesGenerator(seqtype, *inputs)
        tree_generator_obj.build(method=method, threads=threads)
        logging.debug(
            f"Tree building with {phyling.config.avail_tree_methods[method]} was finished in {phyling.config.runtime(func_start)}."
        )
        if len(inputs) == 1:
            tree_generator_obj.get_tree()
        elif top_n_toverr:
            func_start = time.monotonic()
            tree_generator_obj.treeness_over_rcv(threads=threads)
            trees = tree_generator_obj.get_tree(n=top_n_toverr)
            inputs = tree_generator_obj.get_mfa(n=top_n_toverr)
            toverrs = tree_generator_obj.get_toverr(n=top_n_toverr)
            logging.debug(f"Treeness estimation and filtering was finished in {phyling.config.runtime(func_start)}.")
            # Output the name of the selected MSA
            output_selected_trees = output / "top_toverr_trees.tsv"
            with open(output_selected_trees, "w") as f:
                for input, toverr in zip(inputs, toverrs):
                    print(input, toverr, sep="\t", file=f)
            logging.info(f"File name of the selected markers is output to {output_selected_trees}")

    if concat:
        logging.info("Concatenate selected MSAs...")
        func_start = time.monotonic()
        alignmentList = []
        for file in inputs:
            alignment = AlignIO.read(file, format="fasta")
            alignment.annotations["seqtype"] = seqtype
            alignment.annotations["seqname"] = file.name
            alignmentList.append(alignment)
        concat_alignments = concatenate_fasta(samples, alignmentList, threads=threads)
        logging.debug(f"Fasta concatenation was finished in {phyling.config.runtime(func_start)}.")
        # Output the concatenated MSA to file
        concat_file, _ = output_concat_file(output, concat_alignments)
        inputs = [concat_file]
        logging.info("Use the concatednated fasta to generate final tree.")
        func_start = time.monotonic()
        tree_generator_obj = TreesGenerator(seqtype, *inputs)
        tree_generator_obj.build(method=method, threads=threads)
        logging.debug(
            f"Tree building with {phyling.config.avail_tree_methods[method]} was finished in {phyling.config.runtime(func_start)}."
        )
        final_tree = tree_generator_obj.get_tree()
    else:
        func_start = time.monotonic()
        final_tree = run_astral(trees)
        logging.debug(f"Consensus tree estimation with ASTRAL was finished in {phyling.config.runtime(func_start)}.")

    return final_tree


def fill_missing_taxon(taxonList: list, alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
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


def concatenate_fasta(taxonList: list, alignmentList: list[MultipleSeqAlignment], threads: int) -> MultipleSeqAlignment:
    """Concatenate multiple MSA results into a single fasta."""
    concat_alignments = MultipleSeqAlignment([])
    for sample in taxonList:
        concat_alignments.append(SeqRecord(Seq(""), id=sample, description=""))
    concat_alignments.sort()
    with Pool(threads) as pool:
        alignmentList = pool.starmap(fill_missing_taxon, product([[seq.id for seq in concat_alignments]], alignmentList))
    logging.debug("Filling missing taxon done")
    end = 0
    partition_info = []
    for alignment in alignmentList:
        start = end + 1
        end += alignment.get_alignment_length()
        if "seqtype" in alignment.annotations:
            partition_info.append(f'{alignment.annotations["seqtype"]}, {alignment.annotations["seqname"]}={start}-{end}')
        concat_alignments += alignment
    if "seqtype" in alignment.annotations:
        concat_alignments.annotations["partition"] = partition_info
    for seq in concat_alignments:
        seq.description = ""
    logging.debug("Concatenate fasta done")
    return concat_alignments


def output_concat_file(output: Path, concat_alignments: MultipleSeqAlignment) -> tuple[Path, Path]:
    """Output the concatenated fasta and partition file."""
    concat_file = output / f"concat_alignments.{phyling.config.aln_ext}"
    with open(concat_file, "w") as f:
        SeqIO.write(concat_alignments, f, format="fasta")
    if "partition" in concat_alignments.annotations:
        partition_file = output / f"concat_alignments.{phyling.config.partition_ext}"
        with open(partition_file, "w") as f:
            [print(line, file=f) for line in concat_alignments.annotations["partition"]]
        logging.info(f"Concatenated fasta and partition file is output to {concat_file} and {partition_file}")
    else:
        partition_file = None
        logging.info("The concat_alignments object doesn't have partition info. No partition file output.")
    return concat_file, partition_file


def run_astral(trees: list[Phylo.BaseTree.Tree]) -> Phylo.BaseTree.Tree:
    """Run astral to get consensus tree."""
    if not shutil.which("astral"):
        raise phyling.BinaryNotFoundError(
            "Astral not found. "
            "Please build the C++ version from the source following the instruction on https://github.com/chaoszhang/ASTER"
        )
    logging.info("Run ASTRAL to resolve consensus among multiple trees")
    temp = StringIO()
    Phylo.write(trees, temp, "newick")
    temp.seek(0)
    p = subprocess.Popen(
        ["astral", "/dev/stdin"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    stdout, _ = p.communicate(temp.read())
    temp.close()
    return Phylo.read(StringIO(stdout), "newick")


def phylotree(inputs, input_dir, output, method, figure, concat, top_n_toverr, threads, **kwargs):
    """
    Construct a phylogenetic tree based on the results of multiple sequence alignment (MSA).

    By default the consensus tree method will be employed which use a 50% cutoff to represent the majority of all the
    trees. You can use the -c/--concat option to concatenate the MSA and build a single tree instead.

    By default, the UPGMA algorithm is used for tree construction. Users can switch to the Neighbor Joining or FastTree
    by specifying the -m/--method nj/ft.

    Once the tree is built, an ASCII figure representing the tree will be displayed, and a treefile in Newick format
    will be generated as output. Additionally, users can choose to obtain a matplotlib-style figure using the
    -f/--figure option.
    """
    module_start = time.monotonic()

    if input_dir:
        inputs = [file for file in Path(input_dir).glob(f"*.{phyling.config.aln_ext}")]
    else:
        inputs = [Path(file) for file in inputs if file.endswith(phyling.config.aln_ext)]
        input_dir = {file.parent for file in inputs}

        if len(input_dir) > 1:
            logging.error("Inputs do not came from the same folder. Cannot obtain the checkpoint file.")
            logging.error(input_dir)
            sys.exit(1)
        else:
            input_dir = input_dir.pop()

    logging.info(f"Found {len(inputs)} MSA fasta")

    input_dir = Path(input_dir)

    # Get the sample names from checkpoint file under the inputs folder
    checkpoint = input_dir / ".checkpoint.pkl"
    try:
        with open(checkpoint, "rb") as f:
            samples, _ = pickle.load(f)
    except FileNotFoundError:
        logging.error(
            'Checkpoint file ".checkpoint.pkl" not found. Please rerun the align module again without --from_checkpoint.'
        )
    seqtype = tuple(samples.values())[0].seqtype
    logging.info(f"Inputs are {seqtype} sequences.")

    final_tree = mfa_to_finaltree(
        inputs=inputs,
        output=output,
        method=method,
        seqtype=seqtype,
        samples=list(samples.keys()),
        concat=concat,
        top_n_toverr=top_n_toverr,
        threads=threads,
    )

    Phylo.draw_ascii(final_tree)

    output_tree = output / f"{method}_tree.nw"
    logging.info(f"Output tree to {output_tree}")
    with open(output_tree, "w") as f:
        Phylo.write(final_tree, f, "newick")

    if figure:
        fig, ax = plt.subplots(figsize=(20, 12))
        output_fig = output / f"{method}_tree.png"
        Phylo.draw(final_tree, axes=ax)
        fig.savefig(output_fig)

    logging.debug(f"Tree module finished in {phyling.config.runtime(module_start)}.")
