import contextlib
import csv
import logging
import os
import re
import subprocess
import sys
import tempfile
from copy import deepcopy
from functools import partialmethod
from io import BytesIO
from itertools import product
from multiprocessing.dummy import Pool
from pathlib import Path

import pyhmmer
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from clipkit import clipkit
from tqdm import tqdm

# Disable tqdm progress bar implemented in clipkit
tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)


def concat_Bytes_streams(files: list) -> tuple[BytesIO, list]:
    """
    Create a in-memory BytesIO to hold the concatenated fasta files.

    Attributes
    ----------
    files : list
        A list of fasta files.

    Return
    ------
    BytesIO
        A BytesIO object that can use as a regular bytes stream.
    list
        A list that includes a tuple of start and end index of each fasta file.
    """
    concat_stream = BytesIO()
    seq_count = []
    count = 0
    for file in files:
        start_idx = count
        with open(file, "rb") as f:
            for line in f.readlines():
                concat_stream.write(line)
                if line.startswith(b">"):
                    count += 1
            concat_stream.write(b"\n")
        seq_count.append((start_idx, count))
    concat_stream.seek(0)  # Change the stream position to the start of the stream
    return concat_stream, seq_count


def dict_merge(dicts_list: list[dict]) -> dict:
    result = {}
    for d in dicts_list:
        for k, v in d.items():
            result.setdefault(k, set()).add(v)
    return result


class msa_generator:
    def __init__(self, inputs: list[Path]):
        self._inputs = inputs
        concat_stream, self._seq_count = concat_Bytes_streams(inputs)
        seq_file = pyhmmer.easel.SequenceFile(concat_stream, digital=True)
        # Use the concatnated fasta in order to retrieve sequences by index later
        self._sequences = seq_file.read_block()

    def _load_hmms(self) -> dict[pyhmmer.plan7.HMM]:
        hmms = {}
        for hmm_path in list(self._markerset.iterdir()):
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                hmm = hmm_file.read()
            hmms[hmm.name.decode()] = hmm
        logging.info(f"Found {len(hmms)} hmm markers")
        return hmms

    def _get_cutoffs(self) -> dict:
        cutoffs = {}
        try:
            with open(self._markerset.parent / "scores_cutoff", "r") as f:
                for line in csv.reader(f, delimiter="\t"):
                    if line[0].startswith("#"):
                        continue
                    cutoffs[line[0]] = float(line[1])
            logging.info(f"Found {len(cutoffs)} hmm marker cutoffs")
        except FileNotFoundError:
            logging.warning("HMM cutoff file not found")
        return cutoffs

    def _run_hmmsearch(
        self, sample: str, sequence: pyhmmer.easel.DigitalSequence, cutoffs: dict, evalue: float, threads: int = 4
    ) -> dict:
        results = {}
        for hits in pyhmmer.hmmsearch(self._hmms.values(), sequence, cpus=threads):
            cog = hits.query_name.decode()
            for hit in hits:
                if (cutoffs and hit.score < cutoffs[cog]) or (not cutoffs and hit.evalue > evalue):
                    continue
                results.setdefault(cog, hit.name)
                break  # The first hit in hits is the best hit
        logging.info(f"Hmmsearch on {sample} done")
        return results

    def search(self, markerset: Path, evalue: float, threads: int) -> None:
        self._markerset = markerset
        self._hmms = self._load_hmms()
        cutoffs = self._get_cutoffs()
        if len(cutoffs) == len(self._hmms):
            logging.info("Use HMM cutoff file to determine cutoff")
        else:
            logging.info("Use evalue to determine cutoff")
            cutoffs = None

        search_res = []
        self._kh = pyhmmer.easel.KeyHash()
        for idx, sample in enumerate(self._inputs):
            # Select the sequences of each sample
            sub_sequences = self._sequences[self._seq_count[idx][0] : self._seq_count[idx][1]]
            for seq in sub_sequences:
                # Replace description to taxon name
                seq.description = sample.name.encode()
                # Use a KeyHash to store seq.name/index pairs which can be used to retreive
                # ortholog sequences by SequenceObject[kh[seq.name]]
                self._kh.add(seq.name)
            if threads >= 8:
                pass
            else:
                # Single process mode
                logging.debug(f"Run in single process mode with {threads} threads")
                search_res.append(self._run_hmmsearch(sample.name, sub_sequences, cutoffs, evalue, threads))
        # Multi processes mode
        if threads >= 8:
            processes = threads // 4
            logging.debug(f"Run in multiprocesses mode. {processes} jobs with 4 threads for each are run concurrently")
            with Pool(processes) as pool:
                search_res = pool.starmap(
                    self._run_hmmsearch,
                    [
                        (
                            sample.name,
                            self._sequences[self._seq_count[idx][0] : self._seq_count[idx][1]],
                            cutoffs,
                            evalue,
                        )
                        for idx, sample in enumerate(self._inputs)
                    ],
                )
        self.orthologs = dict_merge(search_res)

    @property
    def filter_orthologs(self):
        try:
            self.orthologs = dict(filter(lambda item: len(item[1]) >= 3, self.orthologs.items()))
        except AttributeError:
            logging.error("No orthologs dictionary found. Please make sure the search function was run successfully")
        logging.info(f"Found {len(self.orthologs)} orthologs shared among at least 3 samples")

    def _run_hmmalign(self, hmm: str, hits: set) -> MultipleSeqAlignment:
        # Create an empty SequenceBlock object to store the sequences of the orthologs
        seqs = pyhmmer.easel.DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino())
        for hit in hits:
            seqs.append(self._sequences[self._kh[hit]])

        # HMMalign the ortholog sequences to the corresponding HMM markers
        hmm_profile = self._hmms[hmm]
        MSA = pyhmmer.hmmalign(hmm_profile, seqs, trim=True)

        # Create an empty MultipleSeqAlignment object to store the alignment results
        alignment = MultipleSeqAlignment([])
        for name, aligned_seq, seq_info in zip(MSA.names, MSA.alignment, MSA.sequences):
            alignment.append(
                SeqRecord(
                    Seq(re.sub(r"[ZzBbXx\*\.]", "-", aligned_seq)),
                    id=seq_info.description.decode(),
                    name=name.decode(),
                    description=seq_info.description.decode(),
                )
            )
        return alignment

    def _run_muscle(self, hmm: str, hits: set, output: str) -> MultipleSeqAlignment:
        with open(f"{output}/{hmm}.faa", "wb+") as f:
            for hit in hits:
                seq = self._sequences[self._kh[hit]].copy()
                seq.name = deepcopy(seq.description)
                seq.description = "".encode()
                seq.write(f)

        _ = subprocess.check_call(
            ["muscle", "-align", f"{output}/{hmm}.faa", "-output", f"{output}/{hmm}.aln.faa", "-threads", "1"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )
        alignment = AlignIO.read(f"{output}/{hmm}.aln.faa", "fasta")
        return alignment

    # Fill sequence with "-" for missing samples
    def _fill_missing_taxon(self, taxonList: list, alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
        missing = set(taxonList) - set([seq.id for seq in alignment])
        for sample in missing:
            alignment.append(
                SeqRecord(
                    Seq("-" * alignment.get_alignment_length()),
                    id=sample,
                    description=sample,
                )
            )
        return alignment

    def _run_clipcit(self, hmm: str, alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
        # Use clipkit to trim MSA alignment
        keepD, _ = clipkit.keep_trim_and_log(
            alignment,
            gaps=0.9,
            mode=clipkit.TrimmingMode("gappy"),
            use_log=False,
            outFile=f"{hmm}.faa",
            complement=False,
            char=clipkit.SeqType("aa"),
        )

        clipkit.check_if_all_sites_were_trimmed(keepD)

        seqList = []
        for seq in keepD.keys():
            seqList.append(SeqRecord(Seq(str(keepD[seq])), id=str(seq), description=""))
        alignment = MultipleSeqAlignment(seqList)
        return alignment

    def align(self, output: Path, method: str, non_trim: bool, concat: bool, threads: int) -> None:
        concat_alignments = {sample.name: "" for sample in self._inputs}

        concat_alignments = MultipleSeqAlignment([])
        for sample in self._inputs:
            concat_alignments.append(SeqRecord(Seq(""), id=sample.name, description=""))
        concat_alignments.sort()

        # Parallelize the MSA step
        logging.info(f"Use {method} for MSA")
        logging.info(f"Use {threads} threads to parallelize MSA")
        with Pool(threads) as pool:
            if method == "muscle":
                with tempfile.TemporaryDirectory() as tempdir:
                    logging.debug(f"Create tempdir at: {tempdir}")
                    alignmentList = pool.starmap(
                        self._run_muscle, [(hmm, hits, tempdir) for hmm, hits in self.orthologs.items()]
                    )
            else:
                alignmentList = pool.starmap(self._run_hmmalign, [(hmm, hits) for hmm, hits in self.orthologs.items()])
            logging.info("MSA done")

            alignmentList = pool.starmap(
                self._fill_missing_taxon, product([[seq.id for seq in concat_alignments]], alignmentList)
            )
            logging.info("Filling missing taxon done")

            # Output the alingment fasta without clipkit trimming
            if not non_trim:
                with open(os.devnull, "w") as temp_out, contextlib.redirect_stdout(temp_out):
                    alignmentList = pool.starmap(
                        self._run_clipcit, zip([hmm for hmm in self.orthologs.keys()], alignmentList)
                    )
                logging.info("Clipkit done")

        for hmm, alignment in zip([hmm for hmm in self.orthologs.keys()], alignmentList):
            output_aa = output / f"{hmm}.faa"
            alignment.sort()
            if concat:
                concat_alignments += alignment
            else:
                with open(output_aa, "w") as f:
                    SeqIO.write(alignment, f, format="fasta")

        if concat:
            output_concat = output / "concat_alignments.faa"
            with open(output_concat, "w") as f:
                SeqIO.write(concat_alignments, f, format="fasta")
            logging.info(f"Output concatenated fasta to {output_concat}")
        else:
            logging.info(f"Output individual fasta to folder {output}")


def main(inputs, input_dir, output, markerset, evalue, method, non_trim, concat, threads, **kwargs):
    """
    The align module performs multiple sequence alignment (MSA) on orthologous
    protein sequences that match the hmm markers across samples.

    Initially, Hmmsearch is used to match the samples against a given markerset and
    report the top hit of each sample for each hmm marker, representing "orthologs"
    across all samples. In order to build a tree, minimum of 3 samples should be
    used. If the bitscore cutoff file is present in the hmms folder, it will be used
    as the cutoff. Otherwise, an evalue of 1e-10 will be used as the default cutoff.

    Sequences corresponding to orthologs found in more than 3 samples are extracted
    from each input. These sequences then undergo MSA with hmmalign or muscle. The
    resulting alignment is further trimmed using clipkit by default. You can use the
    -n/--non_trim option to skip the trimming step.

    By default, the alignment results are output separately for each hmm marker. The
    consensus tree method is applying by default to construct a phylogenetic tree.
    You can use the -c/--concat option to concatenate the aligned sequences by sample
    and build a single tree afterward.
    """
    # If args.input_dir is used to instead of args.inputs
    if input_dir:
        inputs = list(Path(input_dir).iterdir())
    else:
        inputs = [Path(sample) for sample in inputs]
    # Check input files, terminate if less than 3 files
    if len(inputs) < 3:
        logging.error("Should have at least 3 input files")
        sys.exit(1)
    output = Path(output)
    output.mkdir(exist_ok=True)
    # Check if output dir is empty
    if any(output.iterdir()):
        logging.warning(f"Output directory {output} is not empty. Aborted")
        sys.exit(1)
    markerset = Path(markerset)

    msa = msa_generator(inputs)
    msa.search(markerset, evalue=evalue, threads=threads)
    msa.filter_orthologs
    msa.align(output, non_trim=non_trim, method=method, concat=concat, threads=threads)
