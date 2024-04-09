"""The PHYling main pipelines."""
import logging
import shutil
import sys
import time
from pathlib import Path

import phyling.config as config
import phyling.exception as exception
import phyling.libphyling as libphyling
import phyling.phylotree as phylotree
import phyling.utils as utils


def align(
    inputs: str | Path | None,
    input_dir: str | Path | None,
    output: str | Path,
    markerset: str | Path,
    evalue: float = 1e-10,
    method: str = "hmmalign",
    non_trim: bool = False,
    threads: int = 1,
    **kwargs,
):
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
    module_start = time.monotonic()

    # If args.input_dir is used to instead of args.inputs
    if input_dir:
        inputs = tuple(Path(input_dir).iterdir())
    elif inputs:
        inputs = tuple(Path(sample) for sample in inputs)
    else:
        logging.error("No input given.")
    # Check input files, terminate if less than 4 files
    if len(inputs) < 4:
        logging.error("Should have at least 4 input files.")
        sys.exit(1)
    logging.info("Loading sequences...")
    inputs: libphyling.SampleList = libphyling.SampleList(inputs)
    logging.info(f"Found {len(inputs)} samples.")

    markerset = Path(markerset)
    if not markerset.exists():
        markerset = Path(config.cfg_dir, config.default_HMM, markerset, "hmms")
    if not markerset.exists():
        logging.error(f"Markerset folder does not exist {markerset} - did you download BUSCO?")
        sys.exit(1)
    logging.info(f"Loading markerset from {markerset}...")
    hmms = libphyling.HMMMarkerSet(markerset, markerset.parent / "scores_cutoff")

    params = {
        "evalue": evalue,
        "method": method,
        "non_trim": non_trim,
        "inputs_checksum": inputs.checksum,
        "hmms_checksum": hmms.checksum,
    }

    precheck_obj = libphyling.OutputPrecheck(output, params=params, samplelist=inputs)
    precheck_obj.precheck()
    orthologs: libphyling.Orthologs
    inputs, orthologs = precheck_obj.data

    if method == "muscle" and not shutil.which("muscle"):
        logging.error(
            'muscle not found. Please install it through "conda install -c bioconda muscle>=5.1" '
            "or build from the source following the instruction on https://github.com/rcedgar/muscle."
        )
        sys.exit(1)

    orthologs = libphyling.search(hmms, inputs, orthologs=orthologs, evalue=evalue, threads=threads)

    try:
        filtered_orthologs = orthologs.filter(min_taxa=4)
    except exception.EmptyWarning:
        logging.error(
            "All orthologs were gone after filtering. Please confirm whether the inputs have sufficient "
            "number of sample or if the correct HMM markers were being used."
        )
        sys.exit(1)

    filtered_orthologs.map(inputs)

    if inputs.seqtype == config.seqtype_cds:
        pep_msa_list, cds_msa_list = libphyling.align(filtered_orthologs, hmms=hmms, method=method, threads=threads)
    else:
        pep_msa_list = libphyling.align(filtered_orthologs, hmms=hmms, method=method, threads=threads)

    if not non_trim:
        msa_list = (
            libphyling.trim(pep_msa_list, cds_msa_list) if inputs.seqtype == config.seqtype_cds else libphyling.trim(pep_msa_list)
        )
    else:
        msa_list = cds_msa_list if inputs.seqtype == config.seqtype_cds else pep_msa_list

    if not msa_list:
        logging.error("Nothing left after trimming. Please disable its through --non_trim and try again.")
        sys.exit(1)

    if inputs.seqtype == config.seqtype_cds:
        ext = config.cds_aln_ext
    else:
        ext = config.prot_aln_ext
    logging.info(f"Output individual fasta to folder {output}...")
    [
        precheck_obj.output_results(msa, precheck_obj.path, ".".join((hmm.decode(), ext)))
        for msa, hmm in zip(msa_list, filtered_orthologs.keys())
    ]
    logging.info("Done.")

    precheck_obj.data = orthologs
    precheck_obj.save_checkpoint()

    logging.debug(f"Align module finished in {utils.runtime(module_start)}.")


def tree(
    inputs: list[str | Path] | None,
    input_dir: str | Path | None,
    output: str | Path,
    method: str = "ft",
    top_n_toverr: int = 50,
    concat: bool = False,
    partition: str | bool = False,
    figure: bool = False,
    threads: int = 1,
    **kwargs,
):
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
    module_start = time.monotonic()

    if input_dir:
        inputs = [Path(file) for file in Path(input_dir).glob(f"*.{config.aln_ext}")]
    else:
        inputs = [Path(file) for file in inputs if str(file).endswith(config.aln_ext)]
        input_dir: set = {file.parent for file in inputs}

        if len(input_dir) > 1:
            logging.error("The inputs aren't in the same folder, which indicates might come from different analysis.")
            sys.exit(1)
        input_dir: Path = input_dir.pop()

    inputs_checksum = utils.get_multifiles_checksum(inputs)
    logging.info(f"Found {len(inputs)} MSA fasta.")
    samples, seqtype = phylotree.determine_samples_and_seqtype(inputs, input_dir)
    logging.info(f"Inputs are {seqtype} sequences.")

    params = {
        "method": method,
        "top_n_toverr": top_n_toverr,
        "concat": concat,
        "partition": partition,
        "inputs_checksum": inputs_checksum,
    }

    precheck_obj = phylotree.OutputPrecheck(output, params=params)
    rerun = precheck_obj.precheck()
    wrapper = precheck_obj.data
    logging.debug(f"rerun = {rerun}")

    if len(inputs) == 1:
        tree, wrapper = phylotree.single_mfa(
            inputs, output=output, method=method, seqtype=seqtype, threads=threads, rerun=rerun, wrapper=wrapper
        )
    else:
        if concat:
            tree, wrapper = phylotree.concat(
                inputs,
                output=output,
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
            tree, wrapper = phylotree.consensus(
                inputs,
                output=output,
                method=method,
                seqtype=seqtype,
                top_n_toverr=top_n_toverr,
                threads=threads,
                rerun=rerun,
                wrapper=wrapper,
            )

    precheck_obj.data = wrapper
    precheck_obj.save_checkpoint()
    precheck_obj.output_results(tree, figure=figure)

    logging.debug(f"Tree module finished in {utils.runtime(module_start)}.")
