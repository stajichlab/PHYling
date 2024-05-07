"""Library of routines for supporting PHYling process."""

from __future__ import annotations

import csv
import hashlib
import logging
import re
import subprocess
import textwrap
from collections import UserDict
from io import BytesIO, StringIO
from multiprocessing import Manager, Pool
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import AnyStr, Iterable, Sequence

import numpy as np
import pyhmmer
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from clipkit.modes import TrimmingMode
from clipkit.msa import MSA
from pyhmmer.easel import DigitalSequenceBlock

import phyling._abc as _abc
import phyling._utils as _utils
import phyling.config as config
import phyling.exception as exception


class SampleSeqs(_abc.SeqFileWrapperABC, Sequence):
    """A wrapper of pyhmmer.easel.DigitalSequenceBlock with convenient implementation to retrieve sequence by name."""

    def __init__(self, file: str | Path) -> None:
        """Initialize the object and perform sequence preprocessing."""
        super().__init__(file)
        self._kh = pyhmmer.easel.KeyHash()
        seqblock: DigitalSequenceBlock = pyhmmer.easel.SequenceFile(self.path, digital=True).read_block()

        if seqblock.alphabet.is_amino():
            self._set_seqtype(config.seqtype_pep)
            logging.debug(f"{self.name} is {self._seqtype} sequences.")
            self._process_pep_seqs(seqblock)

        elif seqblock.alphabet.is_dna():
            self._set_seqtype(config.seqtype_cds)
            logging.debug(f"{self.name} is {self._seqtype} sequences.")
            self._process_cds_seqs(seqblock)

        else:
            raise exception.SeqtypeError(
                f"{self.name} is rna sequences, which is not supported. Please convert them to DNA first."
            )
        self._scanned = False

    def __repr__(self) -> str:
        """Return the string representation."""
        info = f"seqtype={self.seqtype}, scanned={self.is_scanned}"
        return super().__repr__() % info

    def __len__(self) -> int:
        """Return the number of sequences contain in this object."""
        return len(self._cds_seqs) if self.seqtype == config.seqtype_cds else len(self._pep_seqs)

    def __eq__(self, other: SampleSeqs) -> bool:
        """Return true if the two objects have the same name."""
        super().__eq__(other)
        return self.name == other.name

    def __lt__(self, other: SampleSeqs) -> bool:
        """Return true if the sample name of the current object is alphabatically prior to another object."""
        super().__lt__(other)
        return str(self.name) < str(other.name)

    def __hash__(self) -> int:
        """Return the hash value compute by the sample name + seqtype."""
        return hash(self.name + self.seqtype)

    def __getitem__(self, keys: str | tuple[str, bytes]) -> pyhmmer.easel.DigitalSequence | DigitalSequenceBlock:
        """
        Get a pyhmmer DigitalSequence or DigitalSequenceBlock.

        Example:
            SampleSeqs["pep"] -> pyhmmer.easel.DigitalSequenceBlock of entire peptide sequences.
            SampleSeqs["cds"] -> pyhmmer.easel.DigitalSequenceBlock of entire cds sequences if the input file is cds.
            SampleSeqs["pep", b"seqname"] -> pyhmmer.easel.DigitalSequence of the peptide sequence named "seqname".

        Note that the secondary key should be an bytes.
        """
        if not isinstance(keys, tuple):
            key1, key2 = keys, None
        else:
            if len(keys) < 2:
                key1, key2 = keys[0], None
            elif len(keys) == 2:
                key1, key2 = keys
            else:
                raise IndexError(
                    f"Too many keys are given. {self.__class__.__qualname__} accepts at most 2 keys but {len(keys)} were indexed."
                )

        if key1 == config.seqtype_pep:
            seqs = self._pep_seqs
        elif key1 == config.seqtype_cds:
            if self.seqtype == config.seqtype_cds:
                seqs = self._cds_seqs
            else:
                raise AttributeError("Cannot obtain the CDS sequences from a peptide inputs.")
        else:
            raise KeyError(f'First key only accepts "{config.seqtype_pep}" or "{config.seqtype_cds}" but got "{key1}".')

        return seqs[self._kh[key2]] if key2 else seqs

    @property
    def name(self) -> str:
        """Return the sample name without the extension."""
        return re.sub(r"(\.(aa|pep|cds|fna|faa))?\.(fasta|fas|faa|fna|seq|fa)(\.gz)?", "", self.path.name)

    @property
    def checksum(self) -> str:
        """Compute md5 checksum based on the checksums of the sequences."""
        if self.seqtype == config.seqtype_pep:
            seqs = self._pep_seqs
        else:
            seqs = self._cds_seqs
        checksums = 0
        for seq in seqs:
            checksums += seq.checksum()
        return hashlib.md5(str(checksums).encode()).hexdigest()

    @property
    def is_scanned(self):
        """Return whether the peptide sequences have been hmmsearched or not."""
        return self._scanned

    def scanned(self):
        """Mark the sample is hmmsearched."""
        self._scanned = True

    def _process_pep_seqs(self, seqblock: DigitalSequenceBlock) -> None:
        """Process the peptide sequences and add the sequence names (bytes) to KeyHash."""
        self._pep_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), [])

        for seq in seqblock:
            seq.description = self.name.encode()
            self._pep_seqs.append(seq)
            self._kh.add(seq.name)

    def _process_cds_seqs(self, seqblock: DigitalSequenceBlock) -> None:
        """
        Process the cds sequences.

        First check whether the cds fasta contains invalid length which cannot be divided by 3. The passed sequences than
        being translated into peptides and their sequence names (bytes) are added to KeyHash.
        The number of invalid sequence as well their sequence name will be printed to the log.
        """
        self._cds_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.dna(), [])
        self._pep_seqs = DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), [])

        problematic_seqs_name = []
        original_size = len(seqblock)
        while seqblock:
            cds_seq = seqblock.pop(0)
            cds_seq.description = self.name.encode()
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
                f"In the file {self.name}, {problematic_seqs_size}/{original_size} seqs has invalid length. "
                # "The seq names are listed below:"
            )
            # logging.warning(problematic_seqs_name)


class SampleList(_abc.DataListABC[SampleSeqs]):
    """A wrapper that stores all the SampleSeqs for an analysis."""

    def __init__(self, files: Iterable[str | Path | SampleSeqs]) -> None:
        """Initialize the object and perform seqtype and duplicated name checks."""
        super().__init__(files)
        for item in files:
            if isinstance(item, SampleSeqs):
                self.data.append(item)
            elif isinstance(item, (str, Path)):
                self.data.append(SampleSeqs(item))
            else:
                raise TypeError(f"{self.__class__.__qualname__} only accepts list of str/Path/SampleSeqs")

        self._seqtype = self._check_seqtype()
        self._duplicated_name_check()
        self.sort()

    def __repr__(self) -> str:
        """Return the object representation when being called."""
        return super().__repr__() % f"(seqtype={self.seqtype})"

    def __getitem__(self, key: int | slice | str) -> SampleSeqs | SampleList[SampleSeqs]:
        """Get an InputSeqs object from the collection by index or name."""
        if isinstance(key, str):
            if key not in self._names():
                raise KeyError(f"{key}: Sample not found.")
            return self.data[self._names().index(key)]
        return super().__getitem__(key)

    @property
    def seqtype(self):
        """Return whether the sequences are peptide or DNA."""
        return self._seqtype

    @property
    def checksum(self) -> str:
        """Compute md5 checksum based on the checksums of the collected objects."""
        checksum = "".join([data.checksum for data in self.data])
        return hashlib.md5(checksum.encode()).hexdigest()

    def append(self, item: SampleSeqs) -> None:
        """Add new sample to the collection."""
        if not isinstance(item, SampleSeqs):
            raise TypeError(f"Can only append {SampleSeqs.__qualname__} object but got {type(item)}.")
        if item.seqtype != self.seqtype:
            raise exception.SeqtypeError("Item represents different seqtype.")
        if item.name in self._names():
            raise exception.IdenticalKeyError(f"Sample name {item.name} already exists.")
        self.data.append(item)
        self.sort()

    def pop(self, i=0) -> SampleSeqs:
        """Get the data out of the collection by index."""
        if isinstance(i, int):
            return self.data.pop(i)

    def update(self, other: SampleList) -> list:
        """Compare with the previously hmmsearched SampleList and update those need to be rerun."""
        if not isinstance(other, self.__class__):
            raise TypeError(f"Can only update with {self.__class__.__qualname__} object but got {type(other)}.")
        if other.seqtype != self.seqtype:
            raise exception.SeqtypeError("Item represents different seqtype.")
        droplist = []
        for prev_data in other:
            if prev_data.name not in self._names():
                droplist.append(prev_data.name)
            elif prev_data.checksum == self[prev_data.name].checksum and prev_data.is_scanned:
                self[prev_data.name].scanned()
                logging.info(f"Skip hmmsearch for {self[prev_data.name].path.name} since it was already done.")
            else:
                logging.info(f"Sample {self[prev_data.name].path.name} has changes in its fasta file. Redo hmmsearch.")

        return droplist if droplist else None

    def _check_seqtype(self) -> str:
        """Check whether samples have mixed types of sequences."""
        seqtypes = tuple({data.seqtype for data in self.data})
        if len(seqtypes) > 1:
            raise exception.SeqtypeError("Inputs contains more than one seqtypes.")
        return seqtypes[0]

    def _duplicated_name_check(self) -> None:
        """Check whether samples share the same name. Return a dict that use sample names as key."""
        check_dict = {}
        for data in self.data:
            check_dict.setdefault(data.name, []).append(data.path)

        msgs = []
        for sample, files in check_dict.items():
            if len(files) > 1:
                msg = f"'{sample}': {[str(file) for file in files]}"
                msgs.append(msg)

        if msgs:
            msgs = ",".join(msgs)
            raise exception.IdenticalKeyError(f"The following files share the same sample name: [{msgs}]")


class HMMMarkerSet(_abc.DataListABC[pyhmmer.plan7.HMM]):
    """A wrapper that stores the selected HMM markers."""

    def __init__(self, folder: str | Path, cutoff: str | Path | None = None, *, raise_err: bool = True):
        """Initialize the object and perform seqtype and duplicated name checks."""
        if isinstance(folder, (str, Path)):
            folder = Path(folder)
        else:
            raise TypeError('Argument "folder" only accepts of str or Path.')
        files = tuple(file for file in folder.iterdir())
        super().__init__(files)
        for file in files:
            with pyhmmer.plan7.HMMFile(file) as hmm_profile:
                self.data.append(hmm_profile.read())

        logging.info(f"Found {len(self.data)} hmm markers.")

        self._cutoff = self._load_cutoffs(cutoff, cutoff_exit=raise_err) if cutoff else False

    def __repr__(self) -> str:
        """Return the object representation when being called."""
        return super().__repr__() % f"(have_cutoff={self.have_cutoff})"

    def __getitem__(self, key: int | slice | bytes) -> pyhmmer.plan7.HMM | list[pyhmmer.plan7.HMM]:
        """Get an pyhmmer.plan7.HMM object from the collection by index or name (bytes)."""
        if isinstance(key, bytes):
            if key not in self._names():
                raise KeyError(f"{key}: Sample not found.")
            return self.data[self._names().index(key)]
        return self.data[key]

    @property
    def checksum(self) -> str:
        """Compute md5 checksum based on the checksums of the collected HMM profiles."""
        checksums = 0
        for data in self.data:
            checksums += data.checksum
        return hashlib.md5(str(checksums).encode()).hexdigest()

    @property
    def have_cutoff(self) -> bool:
        """Return whether use the trusted cutoff loaded from the file."""
        return self._cutoff

    def _load_cutoffs(self, file: str | Path, *, cutoff_exit: bool) -> bool:
        """Retrieve the The model-specific bit-score cutoffs for each HMM."""
        evalue_hint = " Will use evalue instead."
        try:
            with open(file) as f:
                for line in csv.reader(f, delimiter="\t"):
                    if line[0].startswith("#"):
                        continue
                    cutoff = float(line[1])
                    self[line[0].encode()].cutoffs.trusted = (cutoff, cutoff)
            logging.info("Load model-specific bit-score cutoff for each hmm marker.")
        except FileNotFoundError:
            msg = "HMM cutoff file not found."
            if cutoff_exit:
                raise FileNotFoundError(msg)
            else:
                logging.warning(msg + evalue_hint)
        except KeyError:
            msg = "HMM cutoff file doesn't match the markerset."
            if cutoff_exit:
                raise KeyError(msg)
            else:
                logging.warning(msg + evalue_hint)

        if not all(hmm.cutoffs.trusted for hmm in self.data):
            logging.warning("HMM cutoff file is incomplete." + evalue_hint)
            return False
        logging.info("Use HMM cutoff file to determine cutoff.")
        return True


class Orthologs(UserDict):
    """A dict-like structure that have additional features for SampleList mapping and sequence querying."""

    def __init__(self, data: dict | None = None) -> None:
        """Initialize the object and perform seqtype and duplicated name checks."""
        self.data: dict[AnyStr, set] = {}
        if data is not None:
            self.update(data)
        self._is_mapped = False
        self._seqtype = None

    def __repr__(self) -> str:
        """Return the object representation when being called."""
        if len(repr(self.data)) < 40:
            data = self.data
        else:
            data = [f"{key}: {value}" for key, value in self.items()]
            if len(self) > 20:
                data = data[:5] + ["..."] + data[-5:]
            data = ",\n".join([textwrap.indent(line, prefix=" " * 4) for line in data])
            data = f"{{\n{data}}}"
        seqtype = f", seqtype={self.seqtype}" if self.is_mapped else ""
        return f"{self.__class__.__qualname__}(number={len(self)}, mapped={self.is_mapped}{seqtype}){data}"

    def __eq__(self, other: Orthologs) -> bool:
        """Return true if the data of two objects are the same."""
        if not isinstance(other, self.__class__):
            raise TypeError(f"Can only compare to {self.__class__.__qualname__} object but got {type(other)}.")
        return self.data == other.data

    def __setitem__(self, key: AnyStr, item: Sequence[str, bytes]) -> None:
        """Store the given hits in a set by key."""
        if not isinstance(item, Sequence):
            raise TypeError(f"Only accepts tuple or list but got {type(item)}.")
        if len(item) != 2:
            raise IndexError(f"Only accepts 2 elements (sample_id, seq_id) but got {len(item)}.")

        sample_id, seq_id = item
        if not isinstance(sample_id, str) or not isinstance(seq_id, bytes):
            raise TypeError(f'Only accepts tuple with ("str", "bytes") but got ("{type(sample_id)}", "{type(seq_id)}").')
        self.data.setdefault(key, set()).add((sample_id, seq_id))

    def __getstate__(self):
        """Remove the mapped SampleList object during serialization."""
        state = self.__dict__.copy()
        if "_pep_seqs" in state.keys():
            del state["_pep_seqs"]
        if "_cds_seqs" in state.keys():
            del state["_cds_seqs"]
        state["_is_mapped"] = False
        return state

    @property
    def is_mapped(self):
        """Return whether the object have been mapped to a SampleList."""
        return self._is_mapped

    @property
    def seqtype(self):
        """Compute md5 checksum based on the checksums of the sequences."""
        return self._seqtype

    def update(self, other: dict, /, **kwds):
        """Update the current data."""

        def setitem(hits: set):
            """Add the item to the current collection."""
            for item in hits:
                self[key] = item

        if isinstance(other, dict):
            for key in other.keys():
                setitem(other[key])
        else:
            for key, value in other:
                setitem(value)
        for key, value in kwds.items():
            setitem(value)

    def query(self, key: bytes, seqtype: str = config.seqtype_pep) -> list[pyhmmer.easel.DigitalSequence]:
        """Query the mapped SampleList object and retrieve all the sequences from a given key."""
        if not hasattr(self, "_pep_seqs"):
            raise AttributeError("Need to map to the SampleList object to retrieve the peptide sequences.")
        if seqtype == config.seqtype_pep:
            return self._pep_seqs[key]
        elif seqtype == config.seqtype_cds:
            if not self._cds_seqs:
                raise AttributeError("Orthologs was mapped to a SampleList with peptide sequences where cds is not available.")
            return self._cds_seqs[key]
        else:
            raise KeyError('seqtype only accepts "pep" or "cds"')

    def map(self, sample_list: SampleList, threads: int = 1) -> None:
        """Map to the SampleList object for sequence retrieval."""
        self._pep_seqs, self._cds_seqs = {}, {}
        for key, item in self.data.items():
            pep_seqs, cds_seqs = [], []
            for sample_id, seq_id in item:
                seq = sample_list[sample_id][config.seqtype_pep, seq_id].copy()
                seq.name, seq.description = seq.description, seq.name
                pep_seqs.append(seq)
                if sample_list.seqtype == config.seqtype_cds:
                    seq = sample_list[sample_id][config.seqtype_cds, seq_id].copy()
                    seq.name, seq.description = seq.description, seq.name
                    cds_seqs.append(seq)
            self._pep_seqs[key] = pep_seqs
            if sample_list.seqtype == config.seqtype_cds:
                self._cds_seqs[key] = cds_seqs
        self._is_mapped = True
        self._seqtype = sample_list.seqtype

    def filter(self, min_taxa: int = 0, droplist: Sequence | None = None) -> Orthologs:
        """Filter the orthologs by a given min taxa number or a list contains the samples that no longer exists."""
        ortho_dict = self.data
        if droplist:
            ortho_dict = {key: set(filter(lambda hit: hit[0] not in droplist, value)) for key, value in ortho_dict.items()}
            logging.info(f'Remove hits corresponding to {", ".join([str(sample) for sample in droplist])} from orthologs.')

        if min_taxa and ortho_dict:
            ortho_dict = dict(filter(lambda item: len(item[1]) >= min_taxa, ortho_dict.items()))
            logging.info(f"{len(ortho_dict)} orthologs shared among at least {min_taxa} samples.")

        if not ortho_dict:
            raise exception.EmptyWarning("None of the ortholog left after filtering.")

        return Orthologs(ortho_dict)


class OutputPrecheck(_abc.OutputPrecheckABC):
    """A class that provides features for input/output precheck, checkpoint loading/saving and final MSA output."""

    folder: Path
    ckp: str

    @classmethod
    def setup(cls, *, folder: Path, ckp: str = config.libphyling_checkpoint) -> None:
        """Setup the class variable."""
        super().setup(folder, ckp)

    @classmethod
    def precheck(cls, params: dict, samplelist: SampleList, *, force_rerun: bool = False) -> tuple[SampleList, Orthologs | None]:
        """Check the output folder and determine what orthologs should be removed and what samples should be rerun."""
        if params.keys() != config.libphyling_precheck_params:
            raise KeyError(f"Params should contain keys {config.libphyling_precheck_params}")
        super().precheck(params, (samplelist), force_rerun=force_rerun)
        if not any(cls.folder.iterdir()):
            return samplelist, None
        results = cls._determine_rerun(*(params, samplelist), *cls.load_checkpoint())
        return results

    @classmethod
    def load_checkpoint(cls) -> tuple[dict, SampleList, Orthologs]:
        """
        Load the checkpoint and retrieve the required params/data to determine the rerun status.

        This should be run before precheck.
        """
        return super().load_checkpoint()

    @classmethod
    def save_checkpoint(cls, params: dict, samplelist: SampleList, orthologs: Orthologs) -> None:
        """Save the parameters and data as a checkpoint for rerun."""
        return super().save_checkpoint(params, samplelist, orthologs)

    @classmethod
    @_utils.check_cls_vars("folder")
    def output_results(cls, msa: MultipleSeqAlignment, filename: str) -> None:
        """Output the MSA results in separate files."""
        msa.sort()
        with open(cls.folder / filename, "w") as f:
            SeqIO.write(msa, f, format="fasta")

    @classmethod
    def _determine_rerun(
        cls, cur_params: dict, cur_samples: SampleList, prev_params: dict, prev_samples: SampleList, prev_orthologs: Orthologs
    ) -> tuple[SampleList, Orthologs]:
        """Define the actions that need to do when found a checkpoint file in the given output folder."""
        if prev_samples.seqtype != cur_samples.seqtype:
            raise SystemExit("Seqtype is changed. Aborted.")
        diff_params = {param[0] for param in set(cur_params.items()) ^ set(prev_params.items())}
        if not diff_params:
            raise SystemExit("Files not changed and parameters are identical to the previous run. Aborted.")
        if "markerset" in diff_params:
            raise SystemExit("Markerset is changed. Aborted.")
        if "markerset_cutoff" in diff_params:
            return cur_samples, None
        rm_files = [x for x in cls.folder.iterdir() if x.is_file() if x.name != cls.ckp]
        droplist = cur_samples.update(prev_samples)
        if droplist:
            logging.info("Remove samples that no longer exist from the orthologs collection retrieved from the checkpoint.")
            cur_orthologs = prev_orthologs.filter(droplist=droplist)
        else:
            cur_orthologs = prev_orthologs

        cls._remove_files(files=rm_files)
        return cur_samples, cur_orthologs


@_utils.timing
def search(
    inputs: SampleList,
    markerset: HMMMarkerSet,
    *,
    orthologs: Orthologs | None = None,
    evalue: float = 1e-10,
    threads: int = 1,
) -> Orthologs:
    """
    Search a database of proteins with a set of HMMs to find matching protein sequences.

    The function first loads the hmm profiles and get the cutoff if there has one. Next, it will determine whether to run
    hmmsearch in multiprocess manner. If we have N threads and N < 8, the hmmsearch will not be parallelized and all N
    threads will be used for each round. If N >= 8, there will be N // 4 jobs being launched parallelly and each use 4
    threads.

    After the hmmsearch, all the results will be rearrange to a dictionary with dict["hmm", ["sample"]].
    """
    logging.info("Identify orthologs using hmmsearch...")
    inputs = [input for input in inputs if not input.is_scanned]
    if not orthologs:
        orthologs = Orthologs({})

    if inputs:
        if threads < 8:
            # Single process mode
            logging.debug(f"Sequential mode with {threads} threads.")
            search_res = []
            for input in inputs:
                # Select the sequences of each sample
                search_res.append(_run_hmmsearch(markerset, input, evalue, threads))
        else:
            # Multi processes mode
            threads_per_process = 4
            processes = threads // threads_per_process
            logging.debug(f"Multiprocesses mode: {processes} jobs with 4 threads for each are run concurrently.")
            with ThreadPool(processes) as pool:
                search_res = pool.starmap(
                    _run_hmmsearch,
                    [
                        (
                            markerset,
                            input,
                            evalue,
                            threads_per_process,
                        )
                        for input in inputs
                    ],
                )

        for res in search_res:
            for hmm, sample_id, seq_id in res:
                orthologs[hmm] = (sample_id, seq_id)
        logging.info("Done.")
    return orthologs


@_utils.timing
def align(
    orthologs: Orthologs, markerset: HMMMarkerSet | None = None, *, method: str = "hmmalign", threads: int = 1
) -> tuple[list[MultipleSeqAlignment]] | tuple[list[MultipleSeqAlignment], list[MultipleSeqAlignment]]:
    """
    Align a set of identify orthologous proteins.

    First the sequences that have been identified as ortholog will be aligned through hmmalign or muscle. Next, do DNA to
    peptide back translation if seqtype of the given orthologs object is cds.
    """
    if not orthologs.is_mapped:
        raise AttributeError("Need to map the orthologs to a SampleList first.")
    logging.info(f"Run multiple sequence alignment on peptides using {method}...")
    if method == "hmmalign" and not markerset:
        raise AttributeError(f'Need to specify a {HMMMarkerSet} object to argument "hmms" When using hmmalign.')

    if threads == 1:
        logging.debug("Sequential mode with 1 thread.")
        if method == "muscle":
            pep_msa_list = [_run_muscle(orthologs.query(hmm, config.seqtype_pep)) for hmm in orthologs.keys()]
        else:
            pep_msa_list = [_run_hmmalign(orthologs.query(hmm, config.seqtype_pep), markerset[hmm]) for hmm in orthologs.keys()]
    else:
        logging.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
        manager = Manager()
        with Pool(threads) as pool:
            sharedmem_pep_seqs = manager.list([orthologs.query(hmm, config.seqtype_pep) for hmm in orthologs.keys()])
            if method == "muscle":
                pep_msa_list = pool.map(_run_muscle, sharedmem_pep_seqs)
            else:
                sharedmem_markerset = manager.list([markerset[hmm] for hmm in orthologs.keys()])
                pep_msa_list = pool.starmap(_run_hmmalign, zip(sharedmem_pep_seqs, sharedmem_markerset))
    logging.info("Done.")

    if orthologs.seqtype == config.seqtype_pep:
        return (pep_msa_list,)

    logging.info(f"{config.seqtype_cds} found. Run back translation...")
    if threads == 1:
        logging.debug("Sequential mode with 1 thread.")
        cds_seqs_list = [_to_seqrecord(orthologs.query(hmm, config.seqtype_cds)) for hmm in orthologs.keys()]
        cds_msa_list = [_bp_mrtrans(pep_msa, cds_seq) for pep_msa, cds_seq in zip(pep_msa_list, cds_seqs_list)]
    else:
        logging.debug(f"Multiprocesses mode: {threads} jobs are run concurrently.")
        with ThreadPool(threads) as pool:
            cds_seqs_list = pool.map(_to_seqrecord, [orthologs.query(hmm, config.seqtype_cds) for hmm in orthologs.keys()])
        with Pool(threads) as pool:
            sharedmem_items = manager.list(zip(pep_msa_list, cds_seqs_list))
            cds_msa_list = pool.starmap(_bp_mrtrans, sharedmem_items)
    logging.info("Done.")

    return pep_msa_list, cds_msa_list


@_utils.timing
def trim(
    pep_msa_list: list[MultipleSeqAlignment],
    cds_msa_list: list[MultipleSeqAlignment] | None = None,
    *,
    threads: int = 1,
) -> list[MultipleSeqAlignment]:
    """Trim the uninformative regions."""
    if threads == 1:
        if cds_msa_list:
            cds_msa_trimmed_list = [_trim_gaps(pep_msa, cds_msa) for pep_msa, cds_msa in zip(pep_msa_list, cds_msa_list)]
        else:
            pep_msa_trimmed_list = [_trim_gaps(pep_msa) for pep_msa in pep_msa_list]
    else:
        manager = Manager()
        sharedmem_pep_msa = manager.list(pep_msa_list)
        with Pool(threads) as pool:
            if cds_msa_list:
                sharedmem_cds_msa = manager.list(cds_msa_list)
                cds_msa_trimmed_list = pool.starmap(_trim_gaps, zip(sharedmem_pep_msa, sharedmem_cds_msa))
            else:
                pep_msa_trimmed_list = pool.map(_trim_gaps, sharedmem_pep_msa)
    logging.info("Trimming done.")

    return cds_msa_trimmed_list if cds_msa_list else pep_msa_trimmed_list


def _run_hmmsearch(
    hmms: HMMMarkerSet, input: SampleSeqs, evalue: float = 1e-10, threads: int = 4
) -> list[tuple[str, str, bytes]]:
    """
    Run the hmmsearch process using the pyhmmer library.

    This supports multithreaded running. Empirical testing shows that efficiency drops off after more than
    four (4) CPU threads are used so the default supposed to be the best.
    """
    bit_cutoffs = "trusted" if hmms.have_cutoff else None
    r = []
    for hits in pyhmmer.hmmsearch(hmms, input["pep"], E=evalue, bit_cutoffs=bit_cutoffs, cpus=threads):
        hmm = hits.query_name
        for hit in hits:
            if hit.reported:
                r.append((hmm, input.name, hit.name))
                break  # The first hit in hits is the best hit
    logging.info(f"hmmsearch on {input.path.name} is done.")
    input.scanned()
    return r


def _run_hmmalign(seqs: list[pyhmmer.easel.DigitalSequence], hmm_profile: pyhmmer.plan7.HMM) -> MultipleSeqAlignment:
    """Perform an alignment on a set of protein sequences against a target HMM using pyhmmer."""
    seqs_stream = BytesIO()
    [seq.write(seqs_stream) for seq in seqs]
    seqs_stream.seek(0)
    seqs = pyhmmer.easel.SequenceFile(seqs_stream, digital=True).read_block()
    seqs_stream.close()

    # HMMalign the ortholog sequences to the corresponding HMM markers
    msa = pyhmmer.hmmalign(hmm_profile, seqs)

    temp = BytesIO()
    msa.write(temp, "afa")
    temp.seek(0)
    alignment = AlignIO.read(StringIO(temp.read().decode()), "fasta")
    for seq in alignment:
        seq.seq = _substitute_ambiguous_seq(seq)
    temp.close()
    return alignment


def _run_muscle(seqs: list[pyhmmer.easel.DigitalSequence]) -> MultipleSeqAlignment:
    """Run the multiple sequence alignment tool muscle. This assumes muscle v5 CLI interface and options."""
    seqs_stream = BytesIO()
    [seq.write(seqs_stream) for seq in seqs]
    seqs_stream.seek(0)

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
        seq.seq = _substitute_ambiguous_seq(seq)
    seqs_stream.close()
    return alignment


def _to_seqrecord(seqs: list[pyhmmer.easel.DigitalSequence]) -> list[SeqRecord]:
    """Transform the list of pyhmmer.easel.DigitalSequence to Biopython SeqRecord object."""
    seqs_stream = BytesIO()
    [seq.write(seqs_stream) for seq in seqs]
    seqs_stream.seek(0)

    records = list(SeqIO.parse(StringIO(seqs_stream.read().decode()), "fasta"))
    seqs_stream.close()
    return records


def _bp_mrtrans(pep_msa: MultipleSeqAlignment, cds_rec: list[SeqRecord]) -> MultipleSeqAlignment:
    """Implement a transformer of alignments from protein to mRNA/cDNA coordinates."""
    stop_codons = {"TAA", "TAG", "TGA"}
    codon_size = 3
    gap = "-"
    codon_gap = codon_size * gap

    cds_msa = MultipleSeqAlignment([])
    for pep_seq, cds_seq in zip(pep_msa, cds_rec):
        dna_idx = 0
        dna_align = ""
        for align_idx in range(len(pep_seq)):
            codon = cds_seq[dna_idx : dna_idx + codon_size].seq
            if pep_seq[align_idx] == gap or dna_idx >= len(cds_seq):
                dna_align += codon_gap
                if codon not in stop_codons:
                    continue
            else:
                dna_align += codon
            dna_idx += codon_size
        info = {
            "id": cds_seq.id if hasattr(cds_seq, "id") else pep_seq.id,
            "name": cds_seq.name if hasattr(cds_seq, "name") else pep_seq.name,
            "description": cds_seq.description if hasattr(cds_seq, "description") else pep_seq.description,
        }
        cds_msa.append(SeqRecord(Seq(dna_align), **info))
    return cds_msa


def _trim_gaps(pep_msa: MultipleSeqAlignment, cds_msa: MultipleSeqAlignment = None, gaps: int = 0.9) -> MultipleSeqAlignment:
    """Trim some regions on a MSA if the gappyness ratio (gaps/total length) is higher than the argument "gaps"."""
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


def _substitute_ambiguous_seq(seq: Seq) -> Seq:
    """Substitute ambiguous characters in the sequence and make it all uppercase."""
    return Seq(re.sub(r"[ZzBbXx\*\.]", "-", str(seq.seq))).upper()
