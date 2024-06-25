#!/usr/bin/env python3
"""The PHYling CLI menu and sub-process execution."""
from __future__ import annotations

import argparse
import re
import sys
import textwrap
import time
from pathlib import Path

import phyling._internal._config as _config
from phyling import __author__, __name__, __version__, logger
from phyling.pipeline import align, download, tree


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
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help) is None:
            if action.default not in [argparse.SUPPRESS, None, False]:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help


def _menu_download(
    subparser: argparse._SubParsersAction, parent_parser: argparse.ArgumentParser | None = None
) -> argparse.ArgumentParser:
    """Menu for download module."""
    p_download: argparse.ArgumentParser = subparser.add_parser(
        "download",
        parents=[parent_parser] if parent_parser else [],
        formatter_class=_CustomHelpFormatter,
        help="Download HMM markers",
        description=download.__doc__,
    )
    p_download.add_argument("markerset", metavar='HMM markerset or "list"', help="Name of the HMM markerset")
    p_download.set_defaults(func=download)


def _menu_align(
    subparser: argparse._SubParsersAction, parent_parser: argparse.ArgumentParser | None = None
) -> argparse.ArgumentParser:
    """Menu for align module."""
    p_aln: argparse.ArgumentParser = subparser.add_parser(
        "align",
        parents=[parent_parser] if parent_parser else [],
        formatter_class=_CustomHelpFormatter,
        help="Run multiple sequence alignments against orthologs found among samples",
        description=align.__doc__,
    )
    input_type = p_aln.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        dest="inputs",
        metavar=("file", "files"),
        nargs="+",
        help="Query pepetide/cds fasta or gzipped fasta",
    )
    input_type.add_argument(
        "-I",
        "--input_dir",
        dest="inputs",
        metavar="directory",
        type=Path,
        help="Directory containing query pepetide/cds fasta or gzipped fasta",
    )
    p_aln.add_argument(
        "-o",
        "--output",
        metavar="directory",
        type=Path,
        default=f'phyling-align-{time.strftime("%Y%m%d-%H%M%S", time.gmtime())}',
        help="Output directory of the alignment results (default: phyling-align-[YYYYMMDD-HHMMSS] (UTC timestamp))",
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
        help="Hmmsearch reporting threshold (default: %(default)s, only being used when bitscore cutoff file is not available)",
    )
    p_aln.add_argument(
        "-M",
        "--method",
        choices=_config.avail_align_methods,
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
        default=_config.avail_cpus // 4 * 4,
        help="Threads for hmmsearch and the number of parallelized jobs in MSA step. "
        + "Better be multiple of 4 if using more than 8 threads",
    )
    p_aln.set_defaults(func=align)


def _menu_tree(
    subparser: argparse._SubParsersAction, parent_parser: argparse.ArgumentParser | None = None
) -> argparse.ArgumentParser:
    """Menu for tree module."""
    p_tree: argparse.ArgumentParser = subparser.add_parser(
        "tree",
        parents=[parent_parser] if parent_parser else [],
        formatter_class=_CustomHelpFormatter,
        help="Build a phylogenetic tree based on multiple sequence alignment results",
        description=tree.__doc__,
    )
    input_type = p_tree.add_mutually_exclusive_group(required=True)
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
    p_tree.add_argument(
        "-o",
        "--output",
        metavar="directory",
        type=Path,
        default=f'phyling-tree-{time.strftime("%Y%m%d-%H%M%S", time.gmtime())}',
        help="Output directory of the newick treefile (default: phyling-tree-[YYYYMMDD-HHMMSS] (UTC timestamp))",
    )
    p_tree.add_argument(
        "-M",
        "--method",
        choices=_config.avail_tree_methods.keys(),
        default="ft",
        help="Algorithm used for tree building. (default: %(default)s)\n"
        + "Available options:\n"
        + "\n".join(f"{value}: {key}" for key, value in _config.avail_tree_methods.items()),
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
        default=_config.avail_cpus,
        help="Threads for tree construction",
    )
    p_tree.set_defaults(func=tree)


def _menu() -> argparse.ArgumentParser:
    """Menu for this entry point."""
    parser = argparse.ArgumentParser(
        prog=__name__,
        formatter_class=_CustomHelpFormatter,
        description=main.__doc__,
        epilog=f"Written by {__author__}",
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")

    subparsers = parser.add_subparsers()

    parser.add_argument("-V", "--version", action="version", version=__version__)

    _menu_download(subparsers, parent_parser)
    _menu_align(subparsers, parent_parser)
    _menu_tree(subparsers, parent_parser)

    return parser


def main(args: list[str] | None = None) -> int:
    """
    Package will extract phylogenomic markers and build a phylogenetic tree with these.

    PHYling comprises 3 modules - download, align and tree. The download module can be used to download HMM markerset
    from BUSCO. The align module is the core element of this package which generate multiple sequence alignment among
    the orthologs found across samples. The tree module help to build a phylogenetic tree.
    """
    parser = _menu()

    try:
        args = args or sys.argv[1:]
        if not args:
            parser.print_help(sys.stderr)
            raise SystemExit(0)
        args = parser.parse_args(args)

        args = parser.parse_args()

        if hasattr(args, "threads"):
            if args.threads > _config.avail_cpus:
                args.threads = _config.avail_cpus

        if args.verbose:
            logger.setLevel("DEBUG")
            for handler in logger.handlers:
                handler.setLevel("DEBUG")
            logger.debug("Debug mode enabled.")
            logger.debug(vars(args))

        args.func(**vars(args))

    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1

    except Exception as err:
        logger.error(err)
        return 1

    except SystemExit as err:
        if err.code != 0:
            logger.error(err)
            return 1

    return 0


if __name__ == "__main__":
    main()
