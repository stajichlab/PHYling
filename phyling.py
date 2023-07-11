#!/usr/bin/env python3

import argparse
import logging
import sys
import textwrap
from pathlib import Path

import yaml

from types import SimpleNamespace
# from lib.misc import RecursiveRemoveNone
from download import download
from search import search


def main():

    logging.basicConfig(format="%(asctime)s PHYling %(levelname)s %(message)s", level='INFO')
    logger = logging.getLogger()

    # Create namespace object conf
    args = dict()

    args['script_path'] = Path(__file__).resolve().parent
    args['pep_exts'] = ['.aa.fa', '.aa.fasta', '.faa']
    args['cds_exts'] = ['.cds.fasta', '.fna']
    args['database'] = "https://busco-data.ezlab.org/v5/data"
    args['cfg_dir'] = Path.home() / ".phyling"
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

    parent_parser.add_argument(
        '-p', '--parallel', choices=['serial', 'parallel', 'slurm'], default='serial', help='Parallel the tasks by GNU parallel or slurm')
    parent_parser.add_argument('-f', '--config-file', type=Path, help='config file')
    parent_parser.add_argument('-F', '--force', action='store_true', help='Force to re-run even the outputs exist')
    parent_parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode for debug')
    parent_parser.add_argument('-d', '--dry-run', action='store_true', help='Dry-run mode')

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
        'download', parents=[parent_parser],
        help='Download HMM markers',
        description='Download HMM markerset from BUSCO. Input \"list\" to show urls for all available markersets.')
    p_download.add_argument('markerset', metavar='HMM markerset or \"list\"', help='Name of the HMM markerset')
    p_download.add_argument('-o', '--output', type=Path, default='./HMM', help='Output directory to save HMM markerset')
    p_download.set_defaults(func=download)

    p_search = subparsers.add_parser(
        'search', parents=[parent_parser],
        help='Run HMMsearch of markers against proteomes',
        description=f'Search HMM set against the proteomes. Currently supported peptide fasta extentions: {args.pep_exts}')
    p_search.add_argument('-i', '--input', nargs='+', required=True, help='Query pepetide fasta directory')
    p_search.add_argument('-m', '--markerset', type=Path, required=True, help='Name of the HMM markerset')
    p_search.add_argument('-t', '--threads', type=int, default=1, help='Threads for each task (default=1)')
    p_search.add_argument('-E', '--evalue', type=float, default=1e-10, help='HMM reporting threshold (default=1e-10)')
    p_search.add_argument('-o', '--output', type=Path, required=True, help='Output diretory of the search results')
    p_search.set_defaults(func=search)

    p_aln = subparsers.add_parser(
        'aln', parents=[parent_parser],
        help='Construct individual gene alignments and trim',
        description='Construct unaligned fasta files of protein and cds (if provided) and perform multiple alignments')
    p_aln.add_argument('--pep', type=Path, required=True, help='Query pepetide fasta directory')
    p_aln.add_argument('--cds', type=Path, help='Query cds fasta directory')
    p_aln.add_argument('-m', '--marker', type=Path, required=True,
                       help='HMM marker (markers_3.hmmb) created by PHYling download')
    p_aln.add_argument('-t', '--threads', type=int, default=1, help='Threads for each task (default=1)')
    p_aln.add_argument('-s', '--search-dir', type=Path, required=True, help='HMM Search results directory')
    p_aln.add_argument('-o', '--output', type=Path, required=True, help='Output directory of alignment results')
    # p_aln.set_defaults(func=aln)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.parallel == "slurm" and not args.config_file:
        parser.error("config.yml is required when paralleling by slurm!")

    # Parse args from config file
    # if args.config_file:
    #     with open(args.config_file, 'r') as fh:
    #         conf = yaml.safe_load(fh)
    #     conf = RecursiveRemoveNone(conf)
    #     # Some general args not implemented in argparse but will use in slurm
    #     for other_args in ['general', 'slurm']:
    #         if conf.get(other_args):
    #             for k, v in conf.get(other_args).items():
    #                 setattr(args, k, v)
    #     # Insert the subparser args to sys.argv, which will take by parse_args function below
    #     if conf.get(remaining_args[0]):
    #         for k, v in conf.get(remaining_args[0]).items():
    #             sys.argv.insert(2, str(v))
    #             if k.startswith('--'):
    #                 sys.argv.insert(2, k)
    #     else:
    #         logging.warning(f'No params for {remaining_args[0]} found in config file')

    parser.parse_args(namespace=args)

    if args.verbose:
        logger.setLevel('DEBUG')
        for handler in logger.handlers:
            handler.setLevel('DEBUG')

    args.func(args)


if __name__ == "__main__":
    main()