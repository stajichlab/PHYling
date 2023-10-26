"""Library of routines for supporting PHYling process."""
from __future__ import annotations

import csv
import json
import logging
import re
import shutil
import subprocess
import sys
from collections.abc import Iterator
from io import BytesIO, StringIO
from multiprocessing.dummy import Pool
from pathlib import Path

import numpy as np
import pyhmmer
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from clipkit.modes import TrimmingMode
from clipkit.msa import MSA
from pyhmmer.easel import DigitalSequenceBlock

import phyling.config


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

    clipkit_pep_msa = MSA.from_bio_msa(pep_msa, gap_chars="-")
    clipkit_pep_msa.trim(mode=TrimmingMode.gappy, gap_threshold=gaps)
    pep_msa = clipkit_pep_msa.to_bio_msa()

    if cds_msa:
        if clipkit_pep_msa._site_positions_to_trim.size > 0:
            infoList = [{"id": rec.id, "name": rec.name, "description": rec.description} for rec in cds_msa]
            clipkit_cds_msa = MSA.from_bio_msa(cds_msa)
            pep_trimList_expanded = np.expand_dims(clipkit_pep_msa._site_positions_to_trim, axis=1)
            cds_site_positions_to_trim = (pep_trimList_expanded * np.array([3]) + np.array([0, 1, 2])).flatten()
            clipkit_cds_msa.trim(site_positions_to_trim=cds_site_positions_to_trim)
            cds_msa = clipkit_cds_msa.to_bio_msa()
        bio_msa = cds_msa
    else:
        bio_msa = pep_msa

    for new_rec, rec in zip(bio_msa, infoList):
        new_rec.id, new_rec.name, new_rec.description = rec["id"], rec["name"], rec["description"]

    return bio_msa


class msa_generator:
    """Generate a multiple sequence alignment using hmmer or muscle."""

    def __init__(self, inputs: list[Path]):
        """Initialize the MSA generator object and perform sequence check and preprocessing."""
        self._inputs = sorted(inputs)
        self._fastastrip_re = re.compile(r"(\.(aa|pep|cds|fna|faa))?\.(fasta|fas|faa|fna|seq|fa)(\.gz)?")
        self._inputs_basename_check()
        # Create a dict to retrieve the sequence later
        self._kh = pyhmmer.easel.KeyHash()

        self._seq_count = [self._sequence_preprocess(input=input) for input in self._inputs]

    def search(self, markerset: Path, evalue: float, threads: int) -> None:
        """
        Search a database of proteins with a set of HMMs to find matching protein sequences.

        The function first loads the hmm profiles and get the cutoff if there has one. Next, it will determine whether to run
        hmmsearch in multiprocess manner. If we have N threads and N < 8, the hmmsearch will not be parallelized and all N
        threads will be used for each round. If N >= 8, there will be N // 4 jobs being launched parallelly and each use 4
        threads.

        After the hmmsearch, all the results will be rearrange to a dictionary with dict["hmm", ["sample"]].
        """
        self._markerset = markerset
        self._hmms = self._load_hmms()
        cutoffs = self._get_cutoffs()

        if threads < 8:
            # Single process mode
            logging.debug(f"Run in single process mode with {threads} threads")
            search_res = []
            for input in self._inputs:
                # Select the sequences of each sample
                pep_seqs = self._pep_seqs
                search_res.append(self._run_hmmsearch(input.name, pep_seqs, cutoffs, evalue, threads))
        else:
            # Multi processes mode
            processes = threads // 4
            logging.debug(f"Run in multiprocesses mode. {processes} jobs with 4 threads for each are run concurrently")
            with Pool(processes) as pool:
                search_res = pool.starmap(
                    self._run_hmmsearch,
                    [
                        (
                            input.name,
                            self._pep_seqs[self._seq_count[idx][0] : self._seq_count[idx][1]],
                            cutoffs,
                            evalue,
                        )
                        for idx, input in enumerate(self._inputs)
                    ],
                )

        self.orthologs = {}
        for res in search_res:
            for hmm, seqid in res:
                self.orthologs.setdefault(hmm, set()).add(seqid)

    def filter_orthologs(self) -> None:
        """Filter the found sequence hits from an HMM search. Orthlogs with fewer than 3 hits will be discarded."""
        if not hasattr(self, "orthologs"):
            raise AttributeError(
                "No orthologs dictionary found. Please make sure the search function was run successfully"
            )

        self.orthologs = dict(filter(lambda item: len(item[1]) >= 3, self.orthologs.items()))
        logging.info(f"Found {len(self.orthologs)} orthologs shared among at least 3 samples")

        if not self.orthologs:
            logging.error(
                "All entries are gone after filtering. Please confirm whether the inputs contain an insufficient "
                "number of sequences or if the correct HMM markers are being used."
            )
            sys.exit(1)

    def align(self, output: Path, method: str, non_trim: bool, threads: int) -> None:
        """
        Align a set of identify orthologous proteins.

        First the sequences that have been identified as ortholog will be aligned through hmmalign or muscle. Next, do the dna to
        peptide back translation if self._cds_seqs is found (which means the inputs are cds fasta). Next, run the trimming to
        remove uninformative regions (can be switch off). Finally, output each MSA results in separate files.
        """
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
                ext = phyling.config.cds_aln_ext
            else:
                alignmentList = pep_msa_List
                ext = phyling.config.prot_aln_ext

        logging.info(f"Output individual fasta to folder {output}")
        for hmm, alignment in zip([hmm for hmm in self.orthologs.keys()], alignmentList):
            output_mfa = output / f"{hmm}.{ext}"
            alignment.sort()
            with open(output_mfa, "w") as f:
                SeqIO.write(alignment, f, format="fasta")

        output_metadata = output / ".metadata.json"
        with open(output_metadata, "w") as f:
            json.dump(self._generate_metadata(), f)

    def _sequence_preprocess(self, input: Path) -> tuple(int, int):
        """
        Process the sequence, concatenate all the processed sequences and return the start and end indices.

        When accepting the first input (both self._pep_seqs and self._cds_seqs have not yet created), the function determines the
        input type (peptide or DNA) based on the first input.

        For the following rounds, it checks whether the input have the same data type as determined in the first round.

        All the inputs will be processed by the the following steps:
        1. If the input is a peptide fasta, all the seqs will be appended to the self._pep_seqs. A corresponding entry will be
        added to the self._kh.
        2. If the input is a dna fasta, all the seqs will underwent a valid length check and converted to peptide seq. the dna
        and peptide sequence will be appended to self._pep_seqs and self._cds_seqs, respectively. A corresponding entry will also
        be added to the self._kh.
        3. Finally the function will return a tuple representing the start and end indices of the current file. The indices can
        be used to retrieve the sequences for each sample.
        """
        seqblock = pyhmmer.easel.SequenceFile(input, digital=True).read_block()
        sample = self._fastastrip_re.sub(r"", input.name)

        if not (hasattr(self, "_pep_seqs") or hasattr(self, "_cds_seqs")):
            # If this is the first fasta, do the following check
            if seqblock.alphabet.is_amino():
                logging.info("Inputs are peptide sequences")
                self._pep_seqs = seqblock
                [self._add_to_keyhash(sample, seq) for seq in seqblock]
            elif seqblock.alphabet.is_dna():
                logging.info("Inputs are dna sequences")
                self._cds_translation(input, seqblock)
            else:
                logging.error("Inputs are rna sequences, which are not supported. Please convert them to DNA first")
                sys.exit(1)
            start_idx = 0  # The start idx of the first round is 0
        else:
            # From the second round, check whether the input data type is consistent with the seq data already in record
            # First record the current seq length to be the start idx of the next fasta
            start_idx = len(self._pep_seqs)
            if seqblock.alphabet.is_amino() and not hasattr(self, "_cds_seqs"):
                # cds also create pep_seqs, so need to check whether _cds_seqs is not exist
                for seq in seqblock:
                    self._add_to_keyhash(sample, seq)
                    self._pep_seqs.append(seq)
            elif seqblock.alphabet.is_dna() and hasattr(self, "_cds_seqs"):
                self._cds_translation(input, seqblock)
            else:
                logging.error(
                    "Detect both peptide and DNA fasta in inputs. Mix types of inputs are not allowed. Aborted"
                )
                sys.exit(1)
        return (start_idx, len(self._pep_seqs))

    def _add_to_keyhash(self, sample: str, seq: pyhmmer.easel.DigitalSequence) -> None:
        """Assign the sample name to the pyhmmer.easel.DigitalSequence object and add an entry to keyhash."""
        # Replace description to taxon name
        seq.description = sample.encode()
        # Use a KeyHash to store seq.name/index pairs which can be used to retrieve
        # ortholog sequences by SequenceObject[kh[seq.name]]
        self._kh.add(seq.description + b"-" + seq.name)

    def _cds_translation(self, input: Path, cds_seqblock: DigitalSequenceBlock) -> None:
        """Check whether the cds fasta contains invalid length which cannot be divided by 3.

        Pop the each cds sequence from the seqblock and translate into peptide. The cds and peptide sequences of that is
        successful translated will be appended to the self._pep_seqs and self._cds_seqs respectively and record in self._kh.
        The number of invalid sequence as well their sequence name will be printed to the log.
        """
        if not cds_seqblock.alphabet.is_dna:
            logging.critical("Internal error: the msa_generator._cds_translation only take dna DigitalSequenceBlock.")
        if not hasattr(self, "_cds_seqs"):
            self._cds_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.dna(), [])
        if not hasattr(self, "_pep_seqs"):
            self._pep_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), [])
        sample = self._fastastrip_re.sub(r"", input.name)

        problematic_seqs_name = []
        original_size = len(cds_seqblock)
        while cds_seqblock:
            cds_seq = cds_seqblock.pop(0)
            cds_seq.description = sample.encode()
            try:
                self._pep_seqs.append(cds_seq.translate())
                self._cds_seqs.append(cds_seq)
                self._kh.add(cds_seq.description + b"-" + cds_seq.name)
            except ValueError:
                problematic_seqs_name.append(cds_seq.name.decode())

        if problematic_seqs_name:
            problematic_seqs_size = len(problematic_seqs_name)
            problematic_seqs_name = ", ".join(problematic_seqs_name)
            logging.warning(
                f"In the file {input.name}, {problematic_seqs_size}/{original_size} seqs has invalid length. "
                "The seq names are listed below:"
            )
            logging.warning(problematic_seqs_name)

    def _inputs_basename_check(self) -> None:
        """Check whether inputs share the same basename."""
        check_dict = {}
        [
            check_dict.setdefault(self._fastastrip_re.sub(r"", input.name), []).append(str(input))
            for input in self._inputs
        ]
        if any(len(x) > 1 for x in check_dict.values()):
            logging.error("The following files share the same basename:")
            for x in check_dict.values():
                if len(x) > 1:
                    logging.error(", ".join(x))
            sys.exit(1)

    def _load_hmms(self) -> dict[bytes, pyhmmer.plan7.HMM]:
        """Run the pyhmmer steps for loading HMMs for search or alignment."""
        hmms = {}
        for hmm_path in list(self._markerset.iterdir()):
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                hmm = hmm_file.read()
            hmms[hmm.name.decode()] = hmm
        logging.info(f"Found {len(hmms)} hmm markers")
        return hmms

    def _get_cutoffs(self) -> dict[str, float]:
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

        if len(cutoffs) == len(self._hmms):
            logging.info("Use HMM cutoff file to determine cutoff")
        else:
            logging.info("Use evalue to determine cutoff")
            cutoffs = None

        return cutoffs

    def _run_hmmsearch(
        self, sample: str, sequence: pyhmmer.easel.DigitalSequence, cutoffs: dict, evalue: float, threads: int = 4
    ) -> Iterator[str, bytes]:
        """Run the hmmsearch process using the pyhmmer library.

        This supports multithreaded running. Empirical testing shows that efficiency drops off after more than
        four (4) CPU threads are used so the defaults are best here if applied.
        """
        for hits in pyhmmer.hmmsearch(self._hmms.values(), sequence, cpus=threads):
            hmm = hits.query_name.decode()
            for hit in hits:
                if (cutoffs and hit.score < cutoffs[hmm]) or (not cutoffs and hit.evalue > evalue):
                    continue
                yield hmm, f'{self._fastastrip_re.sub(r"", sample)}-'.encode() + hit.name
                break  # The first hit in hits is the best hit
        logging.info(f"hmmsearch on {sample} done")

    def _get_ortholog_seqs(self, hits: set, seqs: DigitalSequenceBlock) -> BytesIO:
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

    def _generate_metadata(self) -> dict:
        """Generate metadata for phylotree module."""
        metadata = {"sample": {}}
        for input in self._inputs:
            sample = self._fastastrip_re.sub(r"", input.name)
            metadata["sample"].setdefault(sample, sample)

        return metadata


def main(inputs, input_dir, output, markerset, evalue, method, non_trim, threads, **kwargs):
    """
    Perform multiple sequence alignment (MSA) on orthologous sequences that match the hmm markers across samples.

    Initially, HMMsearch is used to match the samples against a given markerset and report the top hit of each sample
    for each hmm marker, representing "orthologs" across all samples. In order to build a tree, minimum of 3 samples
    should be used. If the bitscore cutoff file is present in the hmms folder, it will be used as the cutoff. Otherwise,
    an evalue of 1e-10 will be used as the default cutoff.

    Sequences corresponding to orthologs found in more than 3 samples are extracted from each input. These sequences
    then undergo MSA with hmmalign or muscle. The resulting alignments are further trimmed using clipkit by default.
    You can use the -n/--non_trim option to skip the trimming step. Finally, The alignment results are output
    separately for each hmm marker.
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
    msa.align(output, non_trim=non_trim, method=method, threads=threads)
