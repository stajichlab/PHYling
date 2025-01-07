"""Construct a phylogenetic tree by the selected multiple sequence alignment (MSA) results.

By default the consensus tree method will be employed which use a 50% cutoff to represent the majority of all the trees. You
can use the -c/--concat option to concatenate the MSA and build a single tree instead. Note that enable the -p/--partition
option will also output a partition file that compatible to raxml-ng and IQ-TREE.

For the tree building step, the Fasttree will be used as default algorithm. Users can switch to the raxml-ng or
IQ-TREE by specifying the -m/--method raxml/iqtree.

Once the tree is built, an ASCII figure representing the tree will be displayed, and a treefile in Newick format will be
generated as output. Additionally, users can choose to obtain a matplotlib-style figure using the -f/--figure option.
"""

import argparse
import time
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
from Bio import Phylo

from .. import AVAIL_CPUS, logger
from ..libphyling import TreeMethods, TreeOutputFiles
from ..libphyling._utils import Timer, check_threads
from ..libphyling.tree import MFA2TreeList
from .filter import _input_check


def menu(parser: argparse.ArgumentParser) -> None:
    """Menu for tree module."""
    req_args = parser.add_argument_group("Required arguments")
    input_type = req_args.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        dest="inputs",
        metavar=("file", "files"),
        nargs="+",
        type=Path,
        help="Multiple sequence alignment fasta of the markers",
    )
    input_type.add_argument(
        "-I",
        "--input_dir",
        dest="inputs",
        metavar="directory",
        type=Path,
        help="Directory containing multiple sequence alignment fasta of the markers",
    )
    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument(
        "-o",
        "--output",
        metavar="directory",
        type=Path,
        default="phyling-tree-%s" % time.strftime("%Y%m%d-%H%M%S", time.gmtime()),
        help="Output directory of the newick treefile (default: phyling-tree-[YYYYMMDD-HHMMSS] (UTC timestamp))",
    )
    opt_args.add_argument(
        "-M",
        "--method",
        choices=[m.name.lower() for m in TreeMethods],
        default=TreeMethods.FT.name.lower(),
        help="Algorithm used for tree building. (default: %(default)s)\n"
        + "Available options:\n"
        + "\n".join(f"{m.name.lower()}: {m.method}" for m in TreeMethods),
    )
    opt_args.add_argument(
        "-c",
        "--concat",
        action="store_true",
        help="Concatenated alignment results",
    )
    opt_args.add_argument(
        "-p",
        "--partition",
        action="store_true",
        help="Partitioned analysis by sequence. Only works when --concat enabled.",
    )
    opt_args.add_argument("-f", "--figure", action="store_true", help="Generate a matplotlib tree figure")
    opt_args.add_argument(
        "-t",
        "--threads",
        type=int,
        default=AVAIL_CPUS,
        help="Threads for tree construction",
    )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    parser.set_defaults(func=tree)


@Timer.timer
@check_threads
def tree(
    inputs: str | Path | list[str | Path],
    output: str | Path,
    *,
    method: Literal["ft", "raxml", "iqtree"] = "ft",
    concat: bool = False,
    partition: bool = False,
    figure: bool = False,
    threads: int = 1,
    **kwargs,
) -> None:
    """A pipeline that build phylogenetic tree through either FastTree, raxml-ng or IQ-TREE."""

    inputs = _input_check(inputs)
    mfa2treelist = MFA2TreeList(inputs)
    partition = _validate_partition(partition, method, concat)

    if concat:
        concat_tree = mfa2treelist.concat(output=output, partition=partition, threads=threads)
        tree, tree_building_cmd = concat_tree.build(method, output=output, threads=threads, capture_cmd=True)
    else:
        mfa2treelist.build(method, output=output, threads=threads)
        tree = mfa2treelist.get_consensus_tree()

    """Output the tree in newick format and figure."""
    Phylo.draw_ascii(tree)
    output = Path(output)
    output.mkdir(exist_ok=True)
    output_tree = output / TreeOutputFiles.TREE_NW
    logger.info("Output tree to %s", output_tree)
    with open(output_tree, "w") as f:
        strategy = "concatenate" if concat else "consensus"
        f.write(f"# Final tree is built using {TreeMethods[method.upper()].method} with {strategy} strategy.")
        if concat and partition:
            f.write(" Partition mode is enabled.")
            if tree_building_cmd:
                f.write(f"\n# Tree building cmd: {tree_building_cmd}")
        f.write("\n\n")
        Phylo.write(tree, f, "newick")

    if figure:
        fig, ax = plt.subplots(figsize=(20, 12))
        output_fig = output / TreeOutputFiles.TREE_IMG
        logger.info("Output figure to %s", output_fig)
        Phylo.draw(tree, do_show=False, axes=ax)
        fig.savefig(output_fig)


def _validate_partition(partition: bool, method: str, concat: bool) -> bool:
    """Validate and preprocess partition settings."""
    method = method.upper()
    if not isinstance(partition, bool):
        raise ValueError("Argument partition only accepts boolean value.")
    if partition:
        if method not in (TreeMethods.RAXML.name, TreeMethods.IQTREE.name):
            raise ValueError(f"Partition is not allowed with {TreeMethods[method].method}.")
        if not concat:
            raise ValueError("Partition is not allowed in consensus mode.")

    return partition