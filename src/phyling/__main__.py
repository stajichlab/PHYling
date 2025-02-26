#!/usr/bin/env python3
"""The PHYling CLI menu and sub-process execution."""

from __future__ import annotations

import argparse
import re
import sys
import textwrap
import traceback

from . import AUTHOR, VERSION, logger
from . import __doc__ as _module_doc
from . import __name__ as _module_name
from .pipeline.align import __doc__ as _align_doc
from .pipeline.align import menu as _align_menu
from .pipeline.download import __doc__ as _download_doc
from .pipeline.download import menu as _download_menu
from .pipeline.filter import __doc__ as _filter_doc
from .pipeline.filter import menu as _filter_menu
from .pipeline.tree import __doc__ as _tree_doc
from .pipeline.tree import menu as _tree_menu


class _CustomHelpFormatter(argparse.HelpFormatter):
    """Custom help formatter for argparse to enhance text wrapping and default value display.

    This formatter adjusts the text wrapping for better readability and ensures that the default value of an argument is displayed
    in the help message if not already present.
    """

    def _fill_text(self, text, width, indent):
        """Format the class/function docstring with appropriate wrapping."""
        text = [self._whitespace_matcher.sub(" ", paragraph.strip()) for paragraph in text.split("\n\n") if paragraph.strip()]
        return "\n\n".join([textwrap.fill(line, width) for line in text])

    def _split_lines(self, text, width):
        """Enable multi-line display in argument help message."""
        text = [self._whitespace_matcher.sub(" ", line.strip()) for line in text.split("\n") if line.strip()]
        return [wrapped_line for line in text for wrapped_line in textwrap.wrap(line, width)]

    def _get_help_string(self, action):
        """Allow additional message after default parameter displayed."""
        help = action.help
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help) is None:
            if action.default not in [argparse.SUPPRESS, None, False]:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help


def menu() -> argparse.ArgumentParser:
    """Create and configure the argument parser for the command-line interface.

    This function sets up the main parser, subparsers, and their respective arguments for the CLI. The subparsers include
    `download`, `align`, `filter`, and `tree`, each corresponding to a specific functionality within the module.

    The parser is configured with a custom help formatter and includes global options for displaying help and version information.

    Returns:
        argparse.ArgumentParser: The configured argument parser.

    Subparsers:
        download: Provides options for downloading HMM markers.
        align: Provides options for running multiple sequence alignments against orthologs.
        filter: Provides options for filtering the multiple sequence alignment results.
        tree: Provides options for building a phylogenetic tree based on alignment results.

    Example:
        To display help information for the main parser:
        ```bash
        python script.py -h
        ```

        To access help for a specific subcommand:
        ```bash
        python script.py download -h
        ```
    """
    parser = argparse.ArgumentParser(
        prog=_module_name,
        formatter_class=_CustomHelpFormatter,
        description=_module_doc,
        epilog=f"Written by {AUTHOR}",
        add_help=False,
    )

    parent_parser = argparse.ArgumentParser(add_help=False)

    # Add subparsers
    subparsers = parser.add_subparsers(title="Modules")
    parser_download: argparse.ArgumentParser = subparsers.add_parser(
        "download",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Download HMM markers",
        description=_download_doc,
        add_help=False,
    )
    parser_align: argparse.ArgumentParser = subparsers.add_parser(
        "align",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Run multiple sequence alignments against orthologs found among samples",
        description=_align_doc,
        add_help=False,
    )
    parser_filter: argparse.ArgumentParser = subparsers.add_parser(
        "filter",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Filter the multiple sequence alignment results",
        description=_filter_doc,
        add_help=False,
    )
    parser_tree: argparse.ArgumentParser = subparsers.add_parser(
        "tree",
        parents=[parent_parser],
        formatter_class=_CustomHelpFormatter,
        help="Build a phylogenetic tree based on selected multiple sequence alignment results",
        description=_tree_doc,
        add_help=False,
    )

    # Add arguments
    _download_menu(parser_download)
    _align_menu(parser_align)
    _filter_menu(parser_filter)
    _tree_menu(parser_tree)

    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    return parser


def main(args: list[str] | None = None) -> int:
    """Main function to parse command-line arguments and execute the corresponding functionality.

    This function serves as the entry point for the script. It initializes the argument parser, handles command-line arguments,
    configures logging, and executes the associated function. It also manages exceptions gracefully, providing appropriate logging
    and exit codes.

    Args:
        args (list[str] | None): A list of command-line arguments. If None, defaults to `sys.argv[1:]`.

    Returns:
        int: Exit code. Returns 0 on success, 1 on failure or interruption.

    Raises:
        SystemExit: If no arguments are provided or an error occurs during argument parsing.
    """
    parser = menu()

    try:
        args = args or sys.argv[1:]
        if not args:
            parser.print_help(sys.stderr)
            raise SystemExit(0)
        args = parser.parse_args(args)

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
        logger.debug("%s", traceback.format_exc())
        return 1

    except SystemExit as err:
        if err.code != 0:
            logger.error(err)
            logger.debug("%s", traceback.format_exc())
            return err.code

    return 0


if __name__ == "__main__":
    main()
