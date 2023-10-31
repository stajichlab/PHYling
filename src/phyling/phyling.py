#!/usr/bin/env python3
"""The PHYling CLI menu and sub-process execution."""
from __future__ import annotations

import argparse
import logging
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

import phyling.config
from phyling.download import download
from phyling.libphyling import main as search_align
from phyling.phylotree import phylotree


def parser_submodule(parser, parent_parser) -> None:
    """Parser for command line inputs."""
    subparsers = parser.add_subparsers()

    # Each subparser would finally called the corresponding function
    # (function "download" is called at the end in this case)
    p_download = subparsers.add_parser(
        "download",
        parents=[parent_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Download HMM markers",
        description=textwrap.dedent(download.__doc__),
    )
    p_download.add_argument("markerset", metavar='HMM markerset or "list"', help="Name of the HMM markerset")
    p_download.add_argument(
        "-o",
        "--output",
        type=Path,
        default=f"./{phyling.config.default_HMM}",
        help=f'Output directory to save HMM markerset (default="./{phyling.config.default_HMM}")',
    )
    p_download.set_defaults(func=download)

    p_aln = subparsers.add_parser(
        "align",
        parents=[parent_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Run multiple sequence alignments against orthologs found among samples",
        description=textwrap.dedent(search_align.__doc__),
    )
    input_type = p_aln.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        help="Query pepetide/cds fasta or gzipped fasta",
    )
    input_type.add_argument(
        "-I",
        "--input_dir",
        type=Path,
        help="Directory containing query pepetide/cds fasta or gzipped fasta",
    )
    p_aln.add_argument(
        "-o",
        "--output",
        type=Path,
        default="./align",
        help='Output directory of the alignment results (default="./align")',
    )
    p_aln.add_argument(
        "-m",
        "--markerset",
        type=Path,
        required=True,
        help="Directory of the HMM markerset",
    )
    p_aln.add_argument(
        "-E",
        "--evalue",
        type=float,
        default=1e-10,
        help="Hmmsearch reporting threshold (default=1e-10)",
    )
    p_aln.add_argument(
        "-M",
        "--method",
        choices=["hmmalign", "muscle"],
        default="hmmalign",
        help='Program used for multiple sequence alignment (default="hmmalign")',
    )
    p_aln.add_argument(
        "--non_trim",
        action="store_true",
        help="Report non-trimmed alignment results",
    )
    p_aln.add_argument(
        "--from_checkpoint",
        action="store_true",
        help="Load previous hmmsearch results from .checkpoint.pkl in the output directory",
    )
    p_aln.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Threads for hmmsearch and the number of parallelized jobs in MSA step (default=1)",
    )
    p_aln.set_defaults(func=search_align)

    p_tree = subparsers.add_parser(
        "tree",
        parents=[parent_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Build a phylogenetic tree based on multiple sequence alignment results",
        description=textwrap.dedent(phylotree.__doc__),
    )
    input_type = p_tree.add_mutually_exclusive_group(required=True)
    input_type.add_argument("-i", "--inputs", nargs="+", help="Multiple sequence alignment fasta")
    input_type.add_argument(
        "-I",
        "--input_dir",
        type=Path,
        help="Directory containing multiple sequence alignment fasta",
    )
    p_tree.add_argument(
        "-o",
        "--output",
        type=Path,
        default=".",
        help='Output directory of the newick treefile (default=".")',
    )
    p_tree.add_argument(
        "-M",
        "--method",
        choices=["upgma", "nj", "ft"],
        default="upgma",
        help='Algorithm used for tree building (default="upgma")',
    )
    p_tree.add_argument(
        "-c",
        "--concat",
        action="store_true",
        help="Concatenated alignment results",
    )
    p_tree.add_argument("-f", "--figure", action="store_true", help="Generate a matplotlib tree figure")
    p_tree.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Threads for tree construction (default=1)",
    )
    p_tree.set_defaults(func=phylotree)


def main():
    """
    Package will extract phylogenomic markers and build a phylogenetic tree with these.

    PHYling comprises 3 modules - download, align and tree. The download module can be used to download HMM markerset
    from BUSCO. The align module is the core element of this package which generate multiple sequence alignment among
    the orthologs found across samples. The tree module help to build a phylogenetic tree.
    """
    logging.basicConfig(format=f"%(asctime)s {main.__name__} %(levelname)s %(message)s", level="INFO")
    logger = logging.getLogger()

    # Create namespace object conf
    args = dict()

    args["script_path"] = Path(__file__).resolve().parent
    args["database"] = "https://busco-data.ezlab.org/v5/data"
    args["cfg_dir"] = Path.home() / ".phyling"
    args = SimpleNamespace(**args)

    # Create config folder in $HOME/.phyling
    args.cfg_dir.mkdir(exist_ok=True)

    # Implement shared arguments between sub-menu, reference from
    # https://stackoverflow.com/questions/33645859/how-to-add-common-arguments-to-argparse-subcommands
    # Build parent_parser which contains shared arguments, and do not use it directly
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")

    # The real parser for user
    parser = argparse.ArgumentParser(
        prog=main.__name__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(main.__doc__),
        epilog=textwrap.dedent(main._epilog),
    )

    parser.add_argument("-V", "--version", action="version", version=version("phyling"))

    parser_submodule(parser, parent_parser)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser.parse_args(namespace=args)

    if args.verbose:
        logger.setLevel("DEBUG")
        for handler in logger.handlers:
            handler.setLevel("DEBUG")

    logging.debug(args)
    args.func(**vars(args))


main.__name__ = "PHYling"
main._epilog = """
Written by Jason Stajich (jason.stajich[at]ucr.edu or jasonstajich.phd[at]gmail.com).
Rewritten by Cheng-Hung Tsai (chenghung.tsai[at]email.ucr.edu).

Initially written https://github.com/1KFG/Phylogenomics and https://github.com/stajichlab/phyling.
"""

if __name__ == "__main__":
    main()
