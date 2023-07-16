#!/usr/bin/env python3

import argparse
import logging
import sys
import textwrap
from pathlib import Path

import yaml

from types import SimpleNamespace
from download import download
from libphyling import main as search_align


def main():

    logging.basicConfig(format="%(asctime)s PHYling %(levelname)s %(message)s", level='INFO')
    logger = logging.getLogger()

    # Create namespace object conf
    args = dict()

    args["script_path"] = Path(__file__).resolve().parent
    args["database"] = "https://busco-data.ezlab.org/v5/data"
    args["cfg_dir"] = Path.home() / ".phyling"
    args = SimpleNamespace(**args)

    # Create config folder in $HOME/.phyling
    args.cfg_dir.mkdir(exist_ok=True)

    # Load info to namespace object
    try:
        with open(f"{args.script_path}/info.yml", 'r') as fh:
            info = yaml.safe_load(fh)
        info = SimpleNamespace(**info)
    except:
        logging.error("info.yml not found")
        sys.exit(1)

    # Implement shared arguments between sub-menu, reference from https://stackoverflow.com/questions/33645859/how-to-add-common-arguments-to-argparse-subcommands
    # Build parent_parser which contains shared arguments, and do not use it directly
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument('-v', '--verbose', action='store_true', help="Verbose mode for debug")

    _, remaining_args = parent_parser.parse_known_args(namespace=args)

    # The real parser for user
    parser = argparse.ArgumentParser(prog=info.program,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     usage=textwrap.dedent(info.usage),
                                     description=textwrap.dedent(info.description),
                                     epilog=textwrap.dedent(info.epilog))

    parser.add_argument('-V', '--version', action='version', version=info.version)

    subparsers = parser.add_subparsers()

    # Each subparser would finally called the corresponding function (function "download" is called at the end in this case)
    p_download = subparsers.add_parser(
        "download", parents=[parent_parser],
        help="Download HMM markers",
        description="Download HMM markerset from BUSCO. Input \'list\' to show urls for all available markersets.")
    p_download.add_argument('markerset', metavar="HMM markerset or \'list\'", help="Name of the HMM markerset")
    p_download.add_argument('-o', '--output', type=Path, default="./HMM", help="Output directory to save HMM markerset")
    p_download.set_defaults(func=download)

    p_aln = subparsers.add_parser(
        "align", parents=[parent_parser],
        help="Run multiple sequence alignments against orthologs found among samples",
        description="Get orthologs by hmmsearch. Get multiple sequence alignment results from orthologs.")
    input_type = p_aln.add_mutually_exclusive_group(required=True)
    input_type.add_argument('-i', '--inputs', nargs='+', help="Query pepetide fasta. Should have at least 3 samples")
    input_type.add_argument('-I', '--input_dir', type=Path, help="Directory containing at least 3 query pepetide fasta")
    p_aln.add_argument('-o', '--output', type=Path, required=True, help="Output diretory of the alignment results")
    p_aln.add_argument('-m', '--markerset', type=Path, required=True, help="Name of the HMM markerset")
    p_aln.add_argument('-E', '--evalue', type=float, default=1e-10, help="Hmmsearch reporting threshold (default=1e-10)")
    p_aln.add_argument('-n', '--non_trim', action='store_true', help="Report non-clipkit-trimmed alignment results")
    p_aln.add_argument('-t', '--threads', type=int, default=1, help="Threads for hmmsearch (default=1)")
    p_aln.set_defaults(func=search_align)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser.parse_args(namespace=args)

    if args.verbose:
        logger.setLevel('DEBUG')
        for handler in logger.handlers:
            handler.setLevel('DEBUG')

    logging.debug(args)
    args.func(**vars(args))


if __name__ == "__main__":
    main()