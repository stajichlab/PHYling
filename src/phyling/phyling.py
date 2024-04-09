#!/usr/bin/env python3
"""The PHYling CLI menu and sub-process execution."""
from __future__ import annotations

import argparse
import logging
import sys
import textwrap
import time
from pathlib import Path

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

import phyling.config
from phyling.download import download
from phyling.pipeline import align, tree


class _CustomHelpFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = [self._whitespace_matcher.sub(" ", line).strip() for line in text.split("\n\n") if line != ""]
        return "\n\n".join([textwrap.fill(line, width) for line in text])

    def _split_lines(self, text, width):
        text = [self._whitespace_matcher.sub(" ", line).strip() for line in text.split("\n") if line != ""]
        formatted_text = []
        [formatted_text.extend(textwrap.wrap(line, width)) for line in text]
        # The textwrap module is used only for formatting help.
        # Delay its import for speeding up the common usage of argparse.
        return formatted_text

    def _get_help_string(self, action):
        help = action.help
        if "%(default)" not in action.help:
            if action.default not in [argparse.SUPPRESS, None, False]:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help


def parser_submodule(parser: argparse.ArgumentParser, parent_parser: argparse.ArgumentParser) -> None:
    """Parser for command line inputs."""
    subparsers = parser.add_subparsers()

    # Each subparser would finally called the corresponding function
    # (function "download" is called at the end in this case)
    p_download = subparsers.add_parser(
        "download",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Download HMM markers",
        description=download.__doc__,
    )
    p_download.add_argument("markerset", metavar='HMM markerset or "list"', help="Name of the HMM markerset")
    p_download.set_defaults(func=download)

    p_aln = subparsers.add_parser(
        "align",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Run multiple sequence alignments against orthologs found among samples",
        description=align.__doc__,
    )
    input_type = p_aln.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        metavar=("file", "files"),
        nargs="+",
        help="Query pepetide/cds fasta or gzipped fasta",
    )
    input_type.add_argument(
        "-I",
        "--input_dir",
        metavar="directory",
        type=Path,
        help="Directory containing query pepetide/cds fasta or gzipped fasta",
    )
    p_aln.add_argument(
        "-o",
        "--output",
        metavar="directory",
        type=Path,
        default=f'phyling-align-{time.strftime("%Y%m%d-%H%M%S%z", time.localtime())}',
        help="Output directory of the alignment results (default: %(default)s [current timestamp])",
    )
    p_aln.add_argument(
        "-m",
        "--markerset",
        metavar="directory",
        type=Path,
        required=True,
        help="Directory of the HMM markerset",
    )
    p_aln.add_argument(
        "-E",
        "--evalue",
        metavar="float",
        type=float,
        default=1e-10,
        help="Hmmsearch reporting threshold",
    )
    p_aln.add_argument(
        "-M",
        "--method",
        choices=["hmmalign", "muscle"],
        default="hmmalign",
        help="Program used for multiple sequence alignment",
    )
    p_aln.add_argument(
        "--non_trim",
        action="store_true",
        help="Report non-trimmed alignment results",
    )
    p_aln.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Threads for hmmsearch and the number of parallelized jobs in MSA step. "
        + "Better be multiple of 4 if using more than 8 threads",
    )
    p_aln.set_defaults(func=align)

    p_tree = subparsers.add_parser(
        "tree",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Build a phylogenetic tree based on multiple sequence alignment results",
        description=tree.__doc__,
    )
    input_type = p_tree.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        metavar=("file", "files"),
        nargs="+",
        type=Path,
        help="Multiple sequence alignment fasta of the markers",
    )
    input_type.add_argument(
        "-I",
        "--input_dir",
        metavar="directory",
        type=Path,
        help="Directory containing multiple sequence alignment fasta of the markers",
    )
    p_tree.add_argument(
        "-o",
        "--output",
        metavar="directory",
        type=Path,
        default=f'phyling-tree-{time.strftime("%Y%m%d-%H%M%S%z", time.localtime())}',
        help="Output directory of the newick treefile (default: %(default)s [current timestamp])",
    )
    p_tree.add_argument(
        "-M",
        "--method",
        choices=phyling.config.avail_tree_methods.keys(),
        default="ft",
        help="Algorithm used for tree building. (default: %(default)s)\n"
        + "Available options:\n"
        + "\n".join(f"{value}: {key}" for key, value in phyling.config.avail_tree_methods.items()),
    )
    p_tree.add_argument(
        "-n",
        "--top_n_toverr",
        type=int,
        default=50,
        help="Select the top n markers based on their treeness/RCV for final tree building "
        + "(default: %(default)s. Specify 0 to use all markers)",
    )
    p_tree.add_argument(
        "-c",
        "--concat",
        action="store_true",
        help="Concatenated alignment results",
    )
    p_tree.add_argument(
        "-p",
        "--partition",
        choices=["seq", "codon", "seq+codon"],
        help="Create a partition file by sequence or by codon position when --concat enabled. "
        + '"codon" and "seq+codon" only work when inputs are DNA sequences',
    )
    p_tree.add_argument("-f", "--figure", action="store_true", help="Generate a matplotlib tree figure")
    p_tree.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Threads for tree construction",
    )
    p_tree.set_defaults(func=tree)


def main():
    """
    Package will extract phylogenomic markers and build a phylogenetic tree with these.

    PHYling comprises 3 modules - download, align and tree. The download module can be used to download HMM markerset
    from BUSCO. The align module is the core element of this package which generate multiple sequence alignment among
    the orthologs found across samples. The tree module help to build a phylogenetic tree.
    """
    logging.basicConfig(format=f"%(asctime)s {main.__name__} %(levelname)s %(message)s", level="INFO", force=True)
    logger = logging.getLogger()

    # Create config folder in $HOME/.phyling
    phyling.config.cfg_dir.mkdir(exist_ok=True)

    # Implement shared arguments between sub-menu, reference from
    # https://stackoverflow.com/questions/33645859/how-to-add-common-arguments-to-argparse-subcommands
    # Build parent_parser which contains shared arguments, and do not use it directly
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")

    # The real parser for user
    parser = argparse.ArgumentParser(
        prog=main.__name__,
        formatter_class=_CustomHelpFormatter,
        description=main.__doc__,
        epilog=main._epilog,
    )

    parser.add_argument("-V", "--version", action="version", version=version("phyling"))

    parser_submodule(parser, parent_parser)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

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
