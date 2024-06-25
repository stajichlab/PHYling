"""The PHYling main pipelines."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from urllib.error import URLError

import phyling._internal._config as _config
import phyling._internal._libdownload as _libdownload
import phyling._internal._libphyling as _libphyling
import phyling._internal._libtree as _libtree
import phyling._internal._utils as _utils
import phyling.exception as exception
from phyling import logger


@_utils.timing
def download(markerset, **kwargs) -> None:
    """
    Help to download/update BUSCO v5 markerset to a local folder.

    First it checks whether the metadata file is exist under the config folder ~/.phyling. A missing or outdated file
    will trigger the module to download/update the metadata.

    Passing "list" to markerset argument will list all the available/already downloaded markersets. Passing a valid
    name to the markerset argument will download the markerset to the config folder ~/.phyling/HMM.
    """
    try:
        with _libdownload.BuscoMetadata(_config.database, _config.cfg_dir / _config.metadata) as metadata:
            markerset_list = metadata.online + metadata.local
            if markerset == "list":
                width, _ = shutil.get_terminal_size((80, 24))
                col = width // 40
                markerset_list = [markerset_list[x : x + col] for x in range(0, len(markerset_list), col)]
                col_width = max(len(word) for row in markerset_list for word in row) + 3  # padding

                if metadata.online:
                    msg = "Datasets available online:"
                    _libdownload.wrapper(metadata.online, col=col, col_width=col_width, msg=msg)

                if metadata.local:
                    msg = "Datasets available on local:"
                    _libdownload.wrapper(metadata.local, col=col, col_width=col_width, msg=msg)
            elif markerset in metadata.online:
                metadata.download(markerset)
            else:
                logger.error(f'Markerset {markerset} not available. Please check it again with "list" option.')
                sys.exit(1)
    except URLError as e:
        logger.error(e)
        sys.exit(1)
    except FileExistsError as e:
        logger.info(e)


@_utils.timing
def align(
    inputs: str | Path | list[str | Path],
    output: str | Path,
    *,
    markerset: str | Path,
    evalue: float = 1e-10,
    method: str = "hmmalign",
    non_trim: bool = False,
    threads: int = 1,
    **kwargs,
) -> None:
    """
    Perform multiple sequence alignment (MSA) on orthologous sequences that match the hmm markers across samples.

    Initially, hmmsearch is used to match the samples against a given markerset and report the top hit of each sample
    for each hmm marker, representing "orthologs" across all samples. In order to build a tree, minimum of 4 samples
    should be used. If the bitscore cutoff file is present in the hmms folder, it will be used as the cutoff. Otherwise,
    an evalue of 1e-10 will be used as the default cutoff.

    Sequences corresponding to orthologs found in more than 4 samples are extracted from each input. These sequences
    then undergo MSA with hmmalign or muscle. The resulting alignments are further trimmed using clipkit by default.
    You can use the --non_trim option to skip the trimming step. Finally, The alignment results are output separately
    for each hmm marker.
    """
    if isinstance(inputs, list):
        inputs = tuple(Path(sample) for sample in inputs)
    else:
        inputs = tuple(Path(inputs).iterdir())

    # Check input files, terminate if less than 4 files
    if len(inputs) < 4:
        logger.error("Should have at least 4 input files.")
        sys.exit(1)
    logger.info("Loading sequences...")
    try:
        inputs: _libphyling.SampleList = _libphyling.SampleList(inputs)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    logger.info(f"Found {len(inputs)} samples.")

    markerset = Path(markerset)
    if not markerset.exists():
        markerset = Path(_config.cfg_dir, _config.default_HMM, markerset, "hmms")
    if not markerset.exists():
        logger.error(f"Markerset folder does not exist {markerset} - did you download BUSCO?")
        sys.exit(1)
    logger.info(f"Loading markerset from {markerset}...")
    markerset: _libphyling.HMMMarkerSet = _libphyling.HMMMarkerSet(markerset, markerset.parent / "scores_cutoff", raise_err=False)

    if evalue >= 1:
        logger.error(f"Invalid evalue: {evalue}")
        sys.exit(1)
    if method not in _config.avail_align_methods:
        logger.error(f"Invalid method: {method}")
        sys.exit(1)
    if threads > _config.avail_cpus:
        logger.error(f"The maximum number of threads is {_config.avail_cpus} but got {threads}")
        sys.exit(1)

    params = {
        "inputs": inputs.checksum,
        "markerset": markerset.checksum,
        "markerset_cutoff": "markerset cutoff" if markerset.have_cutoff else evalue,
        "method": method,
        "non_trim": non_trim,
    }

    _libphyling.OutputPrecheck.setup(folder=output)
    try:
        inputs, orthologs = _libphyling.OutputPrecheck.precheck(params, inputs)
    except Exception as e:
        logger.error(e)
        sys.exit(1)

    if method == "muscle" and not shutil.which("muscle"):
        logger.error(
            'muscle not found. Please install it through "conda install -c bioconda muscle>=5.1" '
            "or build from the source following the instruction on https://github.com/rcedgar/muscle."
        )
        sys.exit(1)

    orthologs = _libphyling.search(inputs, markerset, orthologs=orthologs, evalue=evalue, threads=threads)

    try:
        filtered_orthologs = orthologs.filter(min_taxa=4)
    except exception.EmptyWarning:
        logger.error(
            "All orthologs were gone after filtering. Please confirm whether the inputs have sufficient "
            "number of sample or if the correct HMM markers were being used."
        )
        sys.exit(1)

    filtered_orthologs.map(inputs)

    try:
        msa_lists = _libphyling.align(filtered_orthologs, markerset, method=method, threads=threads)
    except Exception as e:
        logger.error(e)
        sys.exit(1)

    if not non_trim:
        fianl_msa_list = _libphyling.trim(*msa_lists, threads=threads)
    else:
        fianl_msa_list = msa_lists[1] if inputs.seqtype == _config.seqtype_cds else msa_lists[0]

    if not fianl_msa_list:
        logger.error("Nothing left after trimming. Please disable its through --non_trim and try again.")
        sys.exit(1)

    if inputs.seqtype == _config.seqtype_cds:
        ext = _config.cds_aln_ext
    else:
        ext = _config.prot_aln_ext
    logger.info(f"Output individual fasta to folder {output}...")

    for msa, hmm in zip(fianl_msa_list, filtered_orthologs.keys()):
        _libphyling.OutputPrecheck.output_results(msa, ".".join((hmm.decode(), ext)))

    logger.info("Done.")

    _libphyling.OutputPrecheck.save_checkpoint(params, inputs, orthologs)


@_utils.timing
def tree(
    inputs: str | Path | list[str | Path],
    output: str | Path,
    *,
    method: str = "ft",
    top_n_toverr: int = 50,
    concat: bool = False,
    partition: str | bool = False,
    figure: bool = False,
    threads: int = 1,
    **kwargs,
) -> None:
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
    if isinstance(inputs, list):
        inputs = tuple(Path(file) for file in inputs)
        input_dir = {file.parent for file in inputs}
        if len(input_dir) > 1:
            logger.error("The inputs aren't in the same folder, which indicates might come from different analysis.")
            sys.exit(1)
        input_dir = input_dir.pop()
    else:
        inputs = Path(inputs)
        if inputs.is_file():
            input_dir = inputs.parent
            inputs = (inputs,)
        else:
            input_dir = inputs
            inputs = tuple(file for file in input_dir.glob(f"*.{_config.aln_ext}"))

    inputs_checksum = _utils.get_multifiles_checksum(inputs)
    logger.info(f"Found {len(inputs)} MSA fasta.")
    samples, seqtype = _libtree.determine_samples_and_seqtype(input_dir)

    if method not in ("raxml", "iqtree") and partition:
        logger.warning("Partition is forced to be disabled since it only works when using raxml and iqtree.")
        partition = None

    if top_n_toverr > len(inputs) or top_n_toverr == 0:
        top_n_toverr = len(inputs)

    params = {
        "method": method,
        "top_n_toverr": top_n_toverr,
        "concat": concat,
        "partition": partition,
        "inputs": inputs_checksum,
    }

    _libtree.OutputPrecheck.setup(folder=output, method=method, concat=concat, partition=partition, figure=figure)
    try:
        rerun, wrapper = _libtree.OutputPrecheck.precheck(params=params)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    logger.debug(f"rerun = {rerun}, wrapper = {wrapper}")

    try:
        if len(inputs) == 1:
            tree, wrapper = _libtree.single_mfa(inputs, output, method=method, seqtype=seqtype, threads=threads)
            params.update(top_n_toverr=1, concat=False, partition=None)
        else:
            if concat:
                tree, wrapper = _libtree.concat(
                    inputs,
                    output,
                    method=method,
                    seqtype=seqtype,
                    samples=samples,
                    top_n_toverr=top_n_toverr,
                    partition=partition,
                    threads=threads,
                    rerun=rerun,
                    wrapper=wrapper,
                )
            else:
                tree, wrapper = _libtree.consensus(
                    inputs,
                    output,
                    method=method,
                    seqtype=seqtype,
                    top_n_toverr=top_n_toverr,
                    threads=threads,
                )
            params.update(partition=_libtree.OutputPrecheck.partition)
    except Exception as e:
        logger.error(e)
        sys.exit(1)

    _libtree.OutputPrecheck.save_checkpoint(params, wrapper)
    _libtree.OutputPrecheck.output_results(tree)
