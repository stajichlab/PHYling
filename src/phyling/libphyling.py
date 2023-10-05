"""Library of routines for supporting PHYling process."""
from __future__ import annotations

import csv
import gzip
import logging
import re
import shutil
import subprocess
import sys
from io import BytesIO, StringIO
from itertools import product
from multiprocessing.dummy import Pool
from pathlib import Path

import numpy as np
import pyhmmer
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import phyling.config


def concat_Bytes_streams(files: list) -> tuple[BytesIO, list]:
    """Create a in-memory BytesIO to hold the concatenated fasta files.

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
        f = open(file, "rb")
        if f.read(2) == b"\x1f\x8b":
            f.close()
            f = gzip.open(file, "rb")
        else:
            f.seek(0)
        for line in f.readlines():
            concat_stream.write(line)
            if line.startswith(b">"):
                count += 1
        concat_stream.write(b"\n")
        f.close()
        seq_count.append((start_idx, count))
    concat_stream.seek(0)  # Change the stream position to the start of the stream
    return concat_stream, seq_count


def dict_merge(dicts_list: list[dict]) -> dict:
    """Merge dictionaries with this helper function."""
    result = {}
    for d in dicts_list:
        for k, v in d.items():
            result.setdefault(k, set()).add(v)
    return result


def bp_mrtrans(pep_msa: MultipleSeqAlignment, cds_seqs: list[SeqRecord]) -> MultipleSeqAlignment:
    """Implement a transformer of alignments from protein to mrna/cdna coordinates."""
    stop_codons = {"TAA", "TAG", "TGA"}
    codon_size = 3
    gap = "-"
    codon_gap = codon_size * gap

    cds_msa = MultipleSeqAlignment([])
    for pep_seq, dna_seq in zip(pep_msa, cds_seqs):
        dna_idx = 0
        dna_align = ""
        for align_idx in range(len(pep_seq)):
            codon = dna_seq[dna_idx : dna_idx + codon_size].seq
            if pep_seq[align_idx] == gap or dna_idx >= len(dna_seq):
                dna_align += codon_gap
                if codon not in stop_codons:
                    continue
            else:
                dna_align += codon
            dna_idx += codon_size
        info = {
            "id": dna_seq.id if hasattr(dna_seq, "id") else pep_seq.id,
            "name": dna_seq.name if hasattr(dna_seq, "name") else pep_seq.name,
            "description": dna_seq.description if hasattr(dna_seq, "description") else pep_seq.description,
        }
        cds_msa.append(SeqRecord(Seq(dna_align), **info))
    return cds_msa


def trim_gaps(
    pep_msa: MultipleSeqAlignment, cds_msa: MultipleSeqAlignment = None, gaps: int = 0.9
) -> MultipleSeqAlignment:
    """Trim the locations on MSA if the gappyness ratio (gaps/total length) is higher than a the arg gaps."""
    if gaps > 1 or gaps < 0:
        raise ValueError('The argument "gaps" should be a float between 0 to 1.')
    infoList = [{"id": rec.id, "name": rec.name, "description": rec.description} for rec in pep_msa]
    msa_array = np.array([list(rec) for rec in pep_msa])
    gappyness = (msa_array == "-").mean(axis=0)
    pep_trimList = np.where(gappyness > gaps)[0]
    if pep_trimList.size > 0:
        msa_array = np.delete(msa_array, pep_trimList, axis=1)

    if cds_msa:
        infoList = [{"id": rec.id, "name": rec.name, "description": rec.description} for rec in cds_msa]
        msa_array = np.array([list(rec) for rec in cds_msa])
        if pep_trimList.size > 0:
            cds_trimList = np.concatenate([np.arange(num * 3, (num + 1) * 3) for num in pep_trimList])
            msa_array = np.delete(msa_array, cds_trimList, axis=1)

    return MultipleSeqAlignment(
        [SeqRecord(Seq("".join(rec)), **info) for rec, info in zip(msa_array.tolist(), infoList)]
    )


class msa_generator:
    """Generate a multiple sequence alignment using hmmer or muscle."""

    def __init__(self, inputs: list[Path]):
        """Initialize the MSA generator object."""
        self._inputs = inputs
        concat_stream, self._seq_count = concat_Bytes_streams(inputs)
        self._sequences = pyhmmer.easel.SequenceFile(concat_stream, digital=True).read_block()
        # Use the concatenated fasta in order to retrieve sequences by index later

        # Create dict for sequence retrieval later
        self._kh = pyhmmer.easel.KeyHash()
        self._fastastrip_re = re.compile(r"(\.(aa|pep|cds|fna|faa))?\.(fasta|fas|faa|fna|seq|fa)(\.gz)?")
        for idx, sample in enumerate(self._inputs):
            # Select the sequences of each sample
            sub_sequences = self._sequences[self._seq_count[idx][0] : self._seq_count[idx][1]]
            for seq in sub_sequences:
                # Replace description to taxon name
                seq.description = self._fastastrip_re.sub(r"", sample.name).encode()
                # Use a KeyHash to store seq.name/index pairs which can be used to retrieve
                # ortholog sequences by SequenceObject[kh[seq.name]]
                self._kh.add(seq.name)
                if len(seq.sequence) % 3 != 0:
                    logging.debug(f"seq {seq.name} is len {len(seq.sequence)} truncating")
                    # create a new sequence object which is truncated
                    seq = pyhmmer.easel.DigitalSequence(seq.alphabet, name=seq.name, description=seq.description,
                                                        sequence=seq.sequence[0: 3 * int(len(seq.sequence) / 3)])

        # Check the inputs are peptide or dna sequences
        if self._sequences.alphabet.is_amino():
            logging.info("Inputs are peptide sequences")
            self._pep_seqs = self._sequences
        elif self._sequences.alphabet.is_dna():
            logging.info("Inputs are dna sequences")
            self._pep_seqs = self._sequences.translate()
            self._cds_seqs = self._sequences
        else:
            logging.error(
                "Inputs are rna sequences, which are not a supported format. Please convert them to DNA first"
            )
            sys.exit()

    def _load_hmms(self) -> dict[pyhmmer.plan7.HMM]:
        """Run the pyhmmer steps for loading HMMs for search or alignment."""
        hmms = {}
        for hmm_path in list(self._markerset.iterdir()):
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                hmm = hmm_file.read()
            hmms[hmm.name.decode()] = hmm
        logging.info(f"Found {len(hmms)} hmm markers")
        return hmms

    def _get_cutoffs(self) -> dict:
        """Retrieve the E-value cutoffs for each specific HMM to use in searching."""
        cutoffs = {}
        try:
            with open(self._markerset.parent / "scores_cutoff") as f:
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
        """Run the hmmsearch process using the pyhmmer library.

        This supports multithreaded running. Empirical testing shows that efficiency drops off after more than
        four (4) CPU threads are used so the defaults are best here if applied.

        """
        results = {}
        for hits in pyhmmer.hmmsearch(self._hmms.values(), sequence, cpus=threads):
            cog = hits.query_name.decode()
            for hit in hits:
                if (cutoffs and hit.score < cutoffs[cog]) or (not cutoffs and hit.evalue > evalue):
                    continue
                results.setdefault(cog, hit.name)
                break  # The first hit in hits is the best hit
        logging.info(f"hmmsearch on {sample} done")
        return results

    def search(self, markerset: Path, evalue: float, threads: int) -> None:
        """Search a database of proteins with a set of HMMs to find matching protein sequences."""
        self._markerset = markerset
        self._hmms = self._load_hmms()
        cutoffs = self._get_cutoffs()
        if len(cutoffs) == len(self._hmms):
            logging.info("Use HMM cutoff file to determine cutoff")
        else:
            logging.info("Use evalue to determine cutoff")
            cutoffs = None

        if threads < 8:
            # Single process mode
            search_res = []
            for idx, sample in enumerate(self._inputs):
                # Select the sequences of each sample
                pep_seqs = self._pep_seqs[self._seq_count[idx][0] : self._seq_count[idx][1]]
                logging.debug(f"Run in single process mode with {threads} threads")
                search_res.append(self._run_hmmsearch(sample.name, pep_seqs, cutoffs, evalue, threads))
        else:
            # Multi processes mode
            processes = threads // 4
            logging.debug(f"Run in multiprocesses mode. {processes} jobs with 4 threads for each are run concurrently")
            with Pool(processes) as pool:
                search_res = pool.starmap(
                    self._run_hmmsearch,
                    [
                        (
                            sample.name,
                            self._pep_seqs[self._seq_count[idx][0] : self._seq_count[idx][1]],
                            cutoffs,
                            evalue,
                        )
                        for idx, sample in enumerate(self._inputs)
                    ],
                )
        self.orthologs = dict_merge(search_res)

    def filter_orthologs(self):
        """Filter the found sequence hits from an HMM search."""
        if hasattr(self, "orthologs"):
            self.orthologs = dict(filter(lambda item: len(item[1]) >= 3, self.orthologs.items()))
            logging.info(f"Found {len(self.orthologs)} orthologs shared among at least 3 samples")
        else:
            raise AttributeError(
                "No orthologs dictionary found. Please make sure the search function was run successfully"
            )

    def _get_ortholog_seqs(self, hits: set, seqs: pyhmmer.easel.DigitalSequenceBlock) -> BytesIO:
        stream = BytesIO()
        for hit in hits:
            pep_seq = seqs[self._kh[hit]].copy()
            pep_seq.name, pep_seq.description = pep_seq.description, pep_seq.name
            pep_seq.write(stream)
        stream.seek(0)
        return stream

    def _run_hmmalign(self, hmm: str, hits: set) -> MultipleSeqAlignment:
        """Perform an alignment of a set of protein sequences against a target HMM using pyhmmer."""
        seqs_stream = self._get_ortholog_seqs(hits, self._pep_seqs)
        seqs = pyhmmer.easel.SequenceFile(seqs_stream, digital=True).read_block()

        # HMMalign the ortholog sequences to the corresponding HMM markers
        hmm_profile = self._hmms[hmm]
        MSA = pyhmmer.hmmalign(hmm_profile, seqs)

        temp = BytesIO()
        MSA.write(temp, "afa")
        temp.seek(0)
        alignment = AlignIO.read(StringIO(temp.read().decode()), "fasta")
        for seq in alignment:
            seq.seq = Seq(re.sub(r"[ZzBbXx\*\.]", "-", str(seq.seq)))
        return alignment

    def _run_muscle(self, hits: set) -> MultipleSeqAlignment:
        """Run the multiple sequence alignment tool muscle. This assumes muscle v5 CLI interface and options."""
        seqs_stream = self._get_ortholog_seqs(hits, self._pep_seqs)

        p = subprocess.Popen(
            [
                "muscle",
                "-align",
                "/dev/stdin",
                "-output",
                "/dev/stdout",
                "-threads",
                "1",
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, _ = p.communicate(seqs_stream.read())
        alignment = AlignIO.read(StringIO(stdout.decode()), "fasta")
        for seq in alignment:
            seq.seq = Seq(re.sub(r"[ZzBbXx\*\.]", "-", str(seq.seq)))
        return alignment

    def _fill_missing_taxon(self, taxonList: list, alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
        """Include empty gap-only strings in alignment for taxa lacking an ortholog."""
        missing = set(taxonList) - {seq.id for seq in alignment}
        for sample in missing:
            alignment.append(
                SeqRecord(
                    Seq("-" * alignment.get_alignment_length()),
                    id=sample,
                    description="",
                )
            )
        return alignment

    def align(self, output: Path, method: str, non_trim: bool, concat: bool, threads: int) -> None:
        """Align a set of identify orthologous proteins."""
        # Parallelize the MSA step
        logging.info(f"Use {method} for peptide MSA")
        logging.info(f"Use {threads} threads to parallelize MSA")
        with Pool(threads) as pool:
            if method == "muscle":
                pep_msa_List = pool.map(self._run_muscle, self.orthologs.values())
            else:
                pep_msa_List = pool.starmap(self._run_hmmalign, self.orthologs.items())
            logging.info("Peptide MSA done")

            if hasattr(self, "_cds_seqs"):
                logging.info("cds found. Processing cds sequences...")
                cds_seqs_stream = pool.starmap(
                    self._get_ortholog_seqs, [(hits, self._cds_seqs) for hits in self.orthologs.values()]
                )
                cds_seqs_List = pool.map(
                    lambda x: [record for record in SeqIO.parse(StringIO(x.read().decode()), "fasta")], cds_seqs_stream
                )
                cds_msa_List = pool.starmap(bp_mrtrans, zip(pep_msa_List, cds_seqs_List))
                logging.info("Back translate complete")

            if not non_trim:
                if hasattr(self, "_cds_seqs"):
                    cds_msa_List = pool.starmap(trim_gaps, zip(pep_msa_List, cds_msa_List))
                else:
                    pep_msa_List = pool.map(trim_gaps, pep_msa_List)
                logging.info("Trimming done")

            if hasattr(self, "_cds_seqs"):
                alignmentList = cds_msa_List
            else:
                alignmentList = pep_msa_List

            if concat:
                concat_alignments = MultipleSeqAlignment([])
                for sample in self._inputs:
                    # strip the filename extension from name so that this resembles sp name
                    seqid = self._fastastrip_re.sub(r"", sample.name)
                    concat_alignments.append(SeqRecord(Seq(""), id=seqid, description=""))
                concat_alignments.sort()
                alignmentList = pool.starmap(
                    self._fill_missing_taxon, product([[seq.id for seq in concat_alignments]], alignmentList)
                )
                logging.info("Filling missing taxon done")

        for hmm, alignment in zip([hmm for hmm in self.orthologs.keys()], alignmentList):
            output_aa = output / f"{hmm}.{phyling.config.protein_ext}"
            alignment.sort()
            if concat:
                concat_alignments += alignment
            else:
                with open(output_aa, "w") as f:
                    SeqIO.write(alignment, f, format="fasta")

        if concat:
            output_concat = output / f"concat_alignments.{phyling.config.prot_aln_ext}"
            with open(output_concat, "w") as f:
                SeqIO.write(concat_alignments, f, format="fasta")
            logging.info(f"Output concatenated fasta to {output_concat}")
        else:
            logging.info(f"Output individual fasta to folder {output}")


def main(inputs, input_dir, output, markerset, evalue, method, non_trim, concat, threads, **kwargs):
    """
    Perform multiple sequence alignment (MSA) on orthologous sequences that match the hmm markers across samples.

    Initially, HMMsearch is used to match the samples against a given markerset and report the top hit of each sample
    for each hmm marker, representing "orthologs" across all samples. In order to build a tree, minimum of 3 samples
    should be used. If the bitscore cutoff file is present in the hmms folder, it will be used as the cutoff. Otherwise,
    an evalue of 1e-10 will be used as the default cutoff.

    Sequences corresponding to orthologs found in more than 3 samples are extracted from each input. These sequences
    then undergo MSA with hmmalign or muscle. The resulting alignment is further trimmed using clipkit by default. You
    can use the -n/--non_trim option to skip the trimming step.

    By default, the alignment results are output separately for each hmm marker. The consensus tree method is applying
    by default to construct a phylogenetic tree. You can use the -c/--concat option to concatenate the aligned
    sequences by sample and build a single tree afterward.
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
    if method == "muscle" and not shutil.which("muscle"):
        logging.error(
            'muscle not found. Please install it through "conda install -c bioconda muscle>=5.1" '
            "or build from the source following the instruction on https://github.com/rcedgar/muscle"
        )
        sys.exit(1)

    if not markerset.exists():
        markerset = Path(phyling.config.default_HMM, markerset, "hmms")
    else:
        markerset = Path(markerset)

    if not markerset.exists():
        logging.error(f"Markerset folder does not exist {markerset} - did you download BUSCO?")
        sys.exit(1)

    msa = msa_generator(inputs)
    msa.search(markerset, evalue=evalue, threads=threads)
    msa.filter_orthologs()
    msa.align(output, non_trim=non_trim, method=method, concat=concat, threads=threads)
