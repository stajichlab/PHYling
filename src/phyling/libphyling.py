"""Library of routines for supporting PHYling process."""
from __future__ import annotations

import csv
import gzip
import hashlib
import logging
import pickle
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


def trim_gaps(pep_msa: MultipleSeqAlignment, cds_msa: MultipleSeqAlignment = None, gaps: int = 0.9) -> MultipleSeqAlignment:
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


class InputSeqs:
    """Provide a class that save metadata and seqs which can be easily retrieve from the msa_generator class."""

    def __init__(self, input: Path):
        """Initialize the InputSeqs object and perform md5 checksum and sequence preprocessing."""
        self._input = Path(input)
        self._get_md5()
        self._kh = pyhmmer.easel.KeyHash()
        self._sequence_preprocess()
        self._scanned = False

    def __repr__(self):
        """Return the str of the Pathlib object when being called."""
        return str(self._input)

    def __eq__(self, other):
        """Set the equal magic method for sorting."""
        return self._input == other._input

    def __lt__(self, other):
        """Set the less-than magic method for sorting."""
        return str(self._input) < str(other._input)

    @property
    def basename(self):
        """Return the file basename."""
        return self._input.name

    @property
    def sample(self):
        """Return the sample name without the extension."""
        return re.sub(r"(\.(aa|pep|cds|fna|faa))?\.(fasta|fas|faa|fna|seq|fa)(\.gz)?", "", self._input.name)

    @property
    def md5(self):
        """Return the md5 checksum to determine whether files have been modified."""
        return self._md5

    @property
    def seqtype(self):
        """Report whether the sequences are peptide or DNA."""
        return self._seqtype

    @property
    def pep_seqs(self):
        """Return the entire peptide sequences."""
        return self._pep_seqs

    @property
    def cds_seqs(self):
        """Return the entire cds sequences if the inputs are cds sequences."""
        if hasattr(self, "_cds_seqs"):
            return self._cds_seqs
        else:
            return None

    @property
    def scanned(self):
        """Report whether the peptide sequences have been hmmsearched or not."""
        return self._scanned

    @scanned.setter
    def scanned(self, value: bool):
        if isinstance(value, bool):
            self._scanned = value
        else:
            raise ValueError("Invalid value. The argument scanned only accept boolean value")

    def get_seq_by_name(self, key, target="pep"):
        """Return the specific sequence by its name."""
        if target == "pep":
            return self.pep_seqs[self._kh[key]]
        elif target == "cds":
            return self.cds_seqs[self._kh[key]]
        else:
            raise AttributeError('Argument target only allow "pep" and "cds".')

    def _get_md5(self):
        """Compute md5 checksum of the file."""
        if self._input.is_file():
            f = open(self._input, "rb")
            if f.read(2) == b"\x1f\x8b":
                f.close()
                f = gzip.open(self._input, "rb")
            else:
                f.seek(0)
            self._md5 = hashlib.md5(f.read()).hexdigest()
            f.close()

    def _sequence_preprocess(self) -> None:
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
        seqblock = pyhmmer.easel.SequenceFile(self._input, digital=True).read_block()

        if seqblock.alphabet.is_amino():
            self._seqtype = "peptide"
            logging.debug(f"{self.basename} is peptide sequences")
            self._pep_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), [])
            for seq in seqblock:
                seq.description = self.sample.encode()
                self._pep_seqs.append(seq)
                self._kh.add(seq.name)

        elif seqblock.alphabet.is_dna():
            self._seqtype = "DNA"
            logging.info(f"{self.basename} is DNA sequences")
            self._cds_translation(seqblock)

        else:
            logging.error(f"{self.basename} is rna sequences, which is not supported. Please convert them to DNA first")
            sys.exit(1)

    def _cds_translation(self, cds_seqblock: DigitalSequenceBlock) -> None:
        """Check whether the cds fasta contains invalid length which cannot be divided by 3.

        Pop the each cds sequence from the seqblock and translate into peptide. The cds and peptide sequences of that is
        successful translated will be appended to the self._pep_seqs and self._cds_seqs respectively and record in self._kh.
        The number of invalid sequence as well their sequence name will be printed to the log.
        """
        self._cds_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.dna(), [])
        self._pep_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), [])

        problematic_seqs_name = []
        original_size = len(cds_seqblock)
        while cds_seqblock:
            cds_seq = cds_seqblock.pop(0)
            cds_seq.description = self.sample.encode()
            try:
                self._pep_seqs.append(cds_seq.translate())
                self._cds_seqs.append(cds_seq)
                self._kh.add(cds_seq.name)
            except ValueError:
                problematic_seqs_name.append(cds_seq.name.decode())

        if problematic_seqs_name:
            problematic_seqs_size = len(problematic_seqs_name)
            problematic_seqs_name = ", ".join(problematic_seqs_name)
            logging.warning(
                f"In the file {self.basename}, {problematic_seqs_size}/{original_size} seqs has invalid length. "
                "The seq names are listed below:"
            )
            logging.warning(problematic_seqs_name)


class msa_generator:
    """Generate a multiple sequence alignment using hmmer or muscle."""

    def __init__(self, inputs: list[Path], checkpoint=None):
        """Initialize the MSA generator object and perform sequence check."""
        inputs = sorted([InputSeqs(input) for input in inputs])
        self._check_inputs_seqtype(inputs)
        self._inputs = self._duplicated_name_check(inputs)

    def load_checkpoint(self, checkpoint: Path):
        """
        Load the checkpoint.pkl as intermediate to prevent redoing hmmsearch when adding new samples.

        This should be run before search.
        """
        try:
            with open(checkpoint, "rb") as f:
                prev_inputs, self.orthologs = pickle.load(f)
        except FileNotFoundError:
            logging.error(
                'Checkpoint file ".checkpoint.pkl" not found. Please rerun the align module again without --from_checkpoint.'
            )

        # Replace the InputSeqs objects with the scanned ones if the md5sum are the same
        for sample in tuple(set(self._inputs).intersection(set(prev_inputs))):
            if self._inputs[sample].md5 != prev_inputs[sample].md5:
                logging.info(f"Sample {self._inputs[sample]} has changes in its fasta file. Redo hmmsearch.")
            else:
                self._inputs[sample] = prev_inputs[sample]
                logging.debug(f"Skip hmmsearch for {self._inputs[sample].basename} since it's already done.")

        # Get the list of samples that no longer exist in the new inputs
        drop_List = [sample for sample in tuple(set(prev_inputs) - set(self._inputs))]

        # Remove genes that no longer exist
        if drop_List:
            for genes in self.orthologs.values():
                to_remove = []
                for seqid in genes:
                    if seqid[0] in drop_List:
                        to_remove.append(seqid)
                for seqid in to_remove:
                    genes.remove(seqid)

            logging.info(
                "The following samples are removed from orthologs collection retrieved from the checkpoint "
                "since they're no longer exist in the current inputs:"
            )
            logging.info(", ".join([str(x) for x in drop_List]))

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
        inputs = [input for input in self._inputs.values() if not input.scanned]

        if threads < 8:
            # Single process mode
            logging.debug(f"Run in single process mode with {threads} threads")
            search_res = []
            for input in inputs:
                # Select the sequences of each sample
                search_res.append(self._run_hmmsearch(input, cutoffs, evalue, threads))
        else:
            # Multi processes mode
            threads_per_process = 4
            processes = threads // threads_per_process
            logging.debug(f"Run in multiprocesses mode. {processes} jobs with 4 threads for each are run concurrently")
            with Pool(processes) as pool:
                search_res = pool.starmap(
                    self._run_hmmsearch,
                    [
                        (
                            input,
                            cutoffs,
                            evalue,
                            threads_per_process,
                        )
                        for input in inputs
                    ],
                )

        if not hasattr(self, "orthologs"):
            self.orthologs = {}
        for res in search_res:
            for hmm, seqid in res:
                self.orthologs.setdefault(hmm, set()).add(seqid)

    def save_checkpoint(self, output: Path) -> None:
        """
        Save the inputs and orthologs to a pickle file to prevent redoing hmmsearch when adding new samples.

        This should be run before filter_orthologs.
        """
        if not hasattr(self, "orthologs"):
            raise AttributeError("No orthologs dictionary found. Please make sure the search function was run successfully")

        checkpoint = output / ".checkpoint.pkl"
        with open(checkpoint, "wb") as f:
            pickle.dump((self._inputs, self.orthologs), f)

    def filter_orthologs(self) -> None:
        """Filter the found sequence hits from an HMM search. Orthlogs with fewer than 3 hits will be discarded."""
        if not hasattr(self, "orthologs"):
            raise AttributeError("No orthologs dictionary found. Please make sure the search function was run successfully")

        self.orthologs = dict(filter(lambda item: len(item[1]) >= 3, self.orthologs.items()))
        logging.info(f"Found {len(self.orthologs)} orthologs shared among at least 3 samples")

        if not self.orthologs:
            logging.error(
                "All entries are gone after filtering. Please confirm whether the inputs contain an insufficient "
                "number of sequences or if the correct HMM markers are being used."
            )
            sys.exit(1)

    def align(self, method: str, non_trim: bool, threads: int) -> None:
        """
        Align a set of identify orthologous proteins.

        First the sequences that have been identified as ortholog will be aligned through hmmalign or muscle. Next, do the dna to
        peptide back translation if self._cds_seqs is found (which means the inputs are cds fasta). Next, run the trimming to
        remove uninformative regions (can be switch off).
        """
        if not hasattr(self, "orthologs"):
            raise AttributeError("No orthologs dictionary found. Please make sure the search function was run successfully")

        # Parallelize the MSA step
        logging.info(f"Use {method} for peptide MSA")
        logging.info(f"Use {threads} threads to parallelize MSA")
        with Pool(threads) as pool:
            if method == "muscle":
                pep_msa_List = pool.map(self._run_muscle, self.orthologs.values())
            else:
                pep_msa_List = pool.starmap(self._run_hmmalign, self.orthologs.items())
            logging.info("Peptide MSA done")

            if self._seqtype == "DNA":
                logging.info("cds found. Processing cds sequences...")
                cds_seqs_stream = pool.starmap(self._get_ortholog_seqs, [(hits, "cds") for hits in self.orthologs.values()])
                cds_seqs_List = pool.map(
                    lambda x: [record for record in SeqIO.parse(StringIO(x.read().decode()), "fasta")], cds_seqs_stream
                )
                cds_msa_List = pool.starmap(bp_mrtrans, zip(pep_msa_List, cds_seqs_List))
                logging.info("Back translate complete")

            if not non_trim:
                if self._seqtype == "DNA":
                    cds_msa_List = pool.starmap(trim_gaps, zip(pep_msa_List, cds_msa_List))
                else:
                    pep_msa_List = pool.map(trim_gaps, pep_msa_List)
                logging.info("Trimming done")

        self.pep_msa_List = pep_msa_List
        if self._seqtype == "DNA":
            self.cds_msa_List = cds_msa_List

    def output_msa(self, output: Path) -> None:
        """Output each MSA results in separate files."""
        if not hasattr(self, "pep_msa_List"):
            raise AttributeError("No pep_msa_list found. Please make sure the align function was run successfully")
        if self._seqtype == "DNA":
            alignmentList = self.cds_msa_List
            ext = phyling.config.cds_aln_ext
        else:
            alignmentList = self.pep_msa_List
            ext = phyling.config.prot_aln_ext

        logging.info(f"Output individual fasta to folder {output}")
        for hmm, alignment in zip([hmm for hmm in self.orthologs.keys()], alignmentList):
            output_mfa = output / f"{hmm}.{ext}"
            alignment.sort()
            with open(output_mfa, "w") as f:
                SeqIO.write(alignment, f, format="fasta")

    def _check_inputs_seqtype(self, inputs: list[InputSeqs]):
        """Check whether inputs have mixed type of sequences."""
        seqtypes = tuple({input.seqtype for input in inputs})
        if len(seqtypes) > 1:
            logging.error("Inputs have mixed type of sequences. Please rerun with -v/--verbose to see the details.")
            sys.exit(1)
        else:
            self._seqtype = seqtypes[0]
            logging.info(f"Inputs are {self._seqtype} sequences.")

    def _duplicated_name_check(self, inputs: list[InputSeqs]) -> dict[str, InputSeqs]:
        """Check whether inputs share the same basename. Return a dict that use basename as key."""
        check_dict = {}
        for input in inputs:
            check_dict.setdefault(input.sample, []).append(input)
        if any(len(lis) > 1 for lis in check_dict.values()):
            logging.error("The following files share the same basename:")
            for lis in check_dict.values():
                if len(lis) > 1:
                    logging.error(", ".join([str(x) for x in lis]))
            sys.exit(1)
        else:
            inputs_dict = {key: value[0] for key, value in check_dict.items()}
        return inputs_dict

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

    def _run_hmmsearch(self, input: InputSeqs, cutoffs: dict, evalue: float, threads: int = 4) -> Iterator[str, bytes]:
        """Run the hmmsearch process using the pyhmmer library.

        This supports multithreaded running. Empirical testing shows that efficiency drops off after more than
        four (4) CPU threads are used so the defaults are best here if applied.
        """
        for hits in pyhmmer.hmmsearch(self._hmms.values(), input.pep_seqs, cpus=threads):
            hmm = hits.query_name.decode()
            for hit in hits:
                if (cutoffs and hit.score < cutoffs[hmm]) or (not cutoffs and hit.evalue > evalue):
                    continue
                yield hmm, (input.sample, hit.name)
                break  # The first hit in hits is the best hit
        logging.info(f"hmmsearch on {input.basename} done")
        input.scanned = True

    def _get_ortholog_seqs(self, hits: set, target="pep") -> BytesIO:
        """Return the set of sequences that matched to a hmm profile."""
        stream = BytesIO()
        for hit in hits:
            input_id, seqid = hit
            pep_seq = self._inputs[input_id].get_seq_by_name(key=seqid, target=target).copy()
            pep_seq.name, pep_seq.description = pep_seq.description, pep_seq.name
            pep_seq.write(stream)
        stream.seek(0)
        return stream

    def _run_hmmalign(self, hmm: str, hits: set[tuple[str, bytes]]) -> MultipleSeqAlignment:
        """Perform an alignment of a set of protein sequences against a target HMM using pyhmmer."""
        seqs_stream = self._get_ortholog_seqs(hits)
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

    def _run_muscle(self, hits: set[tuple[str, bytes]]) -> MultipleSeqAlignment:
        """Run the multiple sequence alignment tool muscle. This assumes muscle v5 CLI interface and options."""
        seqs_stream = self._get_ortholog_seqs(hits)

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


def main(inputs, input_dir, output, markerset, evalue, method, non_trim, from_checkpoint, threads, **kwargs):
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
    old_mfa = [file for file in Path(output).glob(f"*.{phyling.config.aln_ext}")]
    if old_mfa and not from_checkpoint:
        logging.warning(f"Output directory {output} is not empty and from_checkpoint option is not enabled. Aborted")
        sys.exit(1)
    else:
        for file in old_mfa:
            file.unlink()

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
    if from_checkpoint:
        msa.load_checkpoint(output / ".checkpoint.pkl")
    msa.search(markerset, evalue=evalue, threads=threads)
    msa.save_checkpoint(output=output)
    msa.filter_orthologs()
    msa.align(method=method, non_trim=non_trim, threads=threads)
    msa.output_msa(output=output)
