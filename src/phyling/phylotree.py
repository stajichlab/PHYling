"""Phylogenetic tree construction methods."""
from __future__ import annotations

import logging
import pickle
import shutil
import subprocess
import sys
from io import StringIO
from itertools import product
from multiprocessing.dummy import Pool
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import phyling.config


class tree_generator:
    """A phylogentic tree generator used in phyling."""

    def __init__(self, seqtype: str, method: str, *file: Path):
        """Initialize the tree generator object."""
        self._file = file
        self._seqtype = seqtype
        self._method = method
        if self._method == "ft" and not shutil.which("FastTree"):
            logging.error(
                'FastTree not found. Please install it through "conda install -c bioconda fasttree"'
            )
            sys.exit(1)
        if len(self._file) > 1 and not shutil.which("astral"):
            logging.error(
                "Astral not found. "
                "Please build the C++ version from the source following the instruction on https://github.com/chaoszhang/ASTER"
            )
            sys.exit(1)

    def get(self, threads: int) -> Phylo.BaseTree.Tree:
        """Run the phylogeny analysis and get the tree object list."""
        logging.info("Start tree building...")
        if len(self._file) == 1:
            final_tree = self._build(self._file[0])
        else:
            logging.debug(f"Run in multiprocesses mode. {threads} jobs are run concurrently")
            with Pool(threads) as pool:
                trees = pool.map(self._build, self._file)
            final_tree = run_astral(trees)
        return final_tree

    def _build(self, file: Path) -> Phylo.BaseTree.Tree:
        if self._method in ["upgma", "nj"]:
            tree = self._with_phylo_module(file)
        elif self._method == "ft":
            tree = self._with_FastTree(file, self._seqtype)
        logging.debug(f"Tree building on {file.name} is done")
        return tree

    def _with_phylo_module(self, file: Path) -> Phylo.BaseTree.Tree:
        """Run the tree calculation using a simple distance method."""
        MSA = AlignIO.read(file, format="fasta")
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, self._method)
        return constructor.build_tree(MSA)

    def _with_VeryFastTree(self, file: Path, seqtype: str, threads: int) -> Phylo.BaseTree.Tree:
        if seqtype == "DNA":
            cmd = ["VeryFastTree", "-nt", "-gamma", "-threads", str(threads), file]
        else:
            cmd = ["VeryFastTree", "-lg", "-gamma", "-threads", str(threads), file]
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

    def _with_FastTree(self, file: Path, seqtype: str) -> Phylo.BaseTree.Tree:
        if seqtype == "DNA":
            cmd = ["FastTree", "-nt", "-gamma", file]
        else:
            cmd = ["FastTree", "-lg", "-gamma", file]
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
        alignmentList = pool.starmap(
            fill_missing_taxon, product([[seq.id for seq in concat_alignments]], alignmentList)
        )
    logging.debug("Filling missing taxon done")
    for alignment in alignmentList:
        concat_alignments += alignment
    for seq in concat_alignments:
        seq.description = ""
    logging.debug("Concatenate fasta done")
    return concat_alignments


def run_astral(trees: list[Phylo.BaseTree.Tree]) -> Phylo.BaseTree.Tree:
    """Run astral to get consensus tree."""
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


def phylotree(inputs, input_dir, output, method, figure, concat, threads, **kwargs):
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
    method_dict = {"upgma": "UPGMA", "nj": "Neighbor Joining", "ft": "FastTree"}
    logging.info(f"Algorithm choose for tree building: {method_dict[method]}")
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

    output = Path(output)
    output.mkdir(exist_ok=True)

    if concat and len(inputs) > 1:
        # Get the concat_alignments
        logging.info("Generate phylogenetic tree the on concatenated fasta")
        alignmentList = []
        for file in inputs:
            alignmentList.append(AlignIO.read(file, format="fasta"))
        concat_alignments = concatenate_fasta(list(samples.keys()), alignmentList, threads=threads)

        # Output the concatenated MSA to file
        output_concat = output / f"concat_alignments.{phyling.config.aln_ext}"
        with open(output_concat, "w") as f:
            SeqIO.write(concat_alignments, f, format="fasta")
        logging.info(f"Concatenated fasta output to {output_concat}")
        inputs = [output_concat]
    else:
        logging.info("Generate phylogenetic tree on all MSA fasta and conclude a majority consensus tree")

    tree_generator_obj = tree_generator(seqtype, method, *inputs)
    final_tree = tree_generator_obj.get(threads=threads)
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
