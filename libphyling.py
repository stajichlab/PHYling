import contextlib
import logging
import os
import re
import sys
import time
import typing
from functools import partialmethod
from io import BytesIO
from itertools import chain
from pathlib import Path

import pyhmmer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from clipkit import clipkit
from pyhmmer.plan7 import HMM, HMMFile
from tqdm import tqdm


# Disable tqdm progress bar implemented in clipkit
tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

class HMMFiles(typing.ContextManager[typing.Iterable[HMM]]):
    def __init__(self, *files: Path) -> None:
        self.stack = contextlib.ExitStack()
        self.hmmfiles = [self.stack.enter_context(HMMFile(f)) for f in files]

    def __enter__(self) -> typing.Iterable[HMM]:
        return chain.from_iterable(self.hmmfiles)

    def __exit__(self, exc_value: object, exc_type: object, traceback: object) -> None:
        self.stack.close()

def concat_Bytes_streams(files: list) -> (BytesIO, list):
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
        A list that recorded the sequence length of each fasta file.
    """
    concat_stream = BytesIO()
    seq_count = []
    for file in files:
        count = 0
        with open(file, 'rb') as f:
            for line in f.readlines():
                concat_stream.write(line)
                if line.startswith(b'>'):
                    count += 1
            concat_stream.write(b'\n')
        seq_count.append(count)
    concat_stream.seek(0)   # Change the stream position to the start of the stream
    return concat_stream, seq_count

def search(
    inputs: list, markerset: Path, seq_file: pyhmmer.easel.SequenceFile, seq_count: list, evalue: float
    ) -> (dict, pyhmmer.easel.KeyHash):
    all_hmms = list(markerset.iterdir())
    logging.info(f"Found {len(all_hmms)} hmm markers")
    
    orthologs = {}
    kh = pyhmmer.easel.KeyHash()
    for idx, count in enumerate(seq_count):
        sequences = seq_file.read_block(count)
        # Use a KeyHash to store seq.name/index pairs which can be used to retreive
        # ortholog sequences by SequenceObject[kh[seq.name]]
        _ = [kh.add(seq.name) for seq in sequences[:]]
        with HMMFiles(*all_hmms) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, sequences, E=evalue):
                cog = hits.query_name.decode()
                for hit in hits:
                    if hit.included:
                        if cog in orthologs:
                            orthologs[cog].add(hit.name)
                        else:
                            orthologs[cog] = set([hit.name])
                        break   # The first hit in hits is the best hit
        logging.info(f"Hmmsearch on {inputs[idx]} done")

        count = 0
        keys_to_remove = []
        for hmm, hits in orthologs.items():
            if len(hits) < 3:
                continue
            count += 1
    logging.info(f"Found {count} orthologs shared among at least 3 samples")
        
    return orthologs, kh

def align(
    orthologs: dict,
    markerset: Path,
    seq_file: pyhmmer.easel.SequenceFile,
    kh: pyhmmer.easel.KeyHash,
    output: Path,
    non_trim: bool
    ) -> None:
    if non_trim:
        logging.info("Output non-clipkit-trimmed multiple sequence alignment results")
    else:
        logging.debug("Multiple sequence alignment results will be trimmed by clipkit")
    
    # Set the seq_file position to the very beginning
    seq_file.rewind()
    sequences = seq_file.read_block()
    
    with open(os.devnull, "w") as temp_out, contextlib.redirect_stdout(temp_out):
        for hmm, hits in orthologs.items():

            if len(hits) < 3:
                continue

            # Create an empty SequenceBlock object to store the sequences of the orthologs
            seqs = pyhmmer.easel.DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino())
            for hit in hits:
                seqs.append(sequences[kh[hit]])

            # HMMalign the ortholog sequences to the corresponding HMM markers
            with HMMFile(markerset / f"{hmm}.hmm") as hmm_file:
                hmm_profile = hmm_file.read()
                MSA = pyhmmer.hmmalign(hmm_profile, seqs)

            # Create an empty MultipleSeqAlignment object to store the alignment results
            alignment = MultipleSeqAlignment([])
            for name, aligned_seq, seq_info in zip(MSA.names, MSA.alignment, MSA.sequences):
                alignment.append(
                    SeqRecord(Seq(re.sub(r"[ZzBbXx\*\.]", "-", aligned_seq)),
                    id=name.decode(),
                    description=seq_info.description.decode())
                    )

            output_aa = output / f"{hmm}.faa"
            # Output the MSA fasta without clipkit trimming
            if non_trim:
                with open(output_aa, 'w') as f:
                    SeqIO.write(MSA, f, format="fasta")
            else:
                # Use clipkit to trim MSA alignment
                clipkit_start_time = time.time()
                
                keepD, trimD = clipkit.keep_trim_and_log(
                    alignment,
                    gaps=0.9,
                    mode=clipkit.TrimmingMode("gappy"),
                    use_log=False,
                    outFile=output_aa,
                    complement=False,
                    char=clipkit.SeqType("aa")
                    )
                
                clipkit.check_if_all_sites_were_trimmed(keepD)
                clipkit.check_if_entry_contains_only_gaps(keepD)
                clipkit.write_keepD(keepD, output_aa, clipkit.FileFormat("fasta"))

    seq_file.close()

def main(inputs, input_dir, output, markerset, evalue, non_trim, **kwargs):
    # If args.input_dir is used to instead of args.inputs
    if input_dir:
        inputs = list(Path(input_dir).iterdir())
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

    # Concatnate all fasta together (in order to retrieve sequences by index later)
    concat_fasta, seq_count = concat_Bytes_streams(inputs)
    seq_file = pyhmmer.easel.SequenceFile(concat_fasta, digital=True)
    orthologs, kh = search(inputs, markerset, seq_file, seq_count, evalue)
    align(orthologs, markerset, seq_file, kh, output, non_trim)
