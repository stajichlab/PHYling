import contextlib
import logging
import typing
from itertools import chain
from pathlib import Path

import pyhmmer
from pyhmmer.plan7 import HMM, HMMFile


class HMMFiles(typing.ContextManager[typing.Iterable[HMM]]):
    def __init__(self, *files: Path) -> None:
        self.stack = contextlib.ExitStack()
        self.hmmfiles = [self.stack.enter_context(HMMFile(f)) for f in files]

    def __enter__(self) -> typing.Iterable[HMM]:
        return chain.from_iterable(self.hmmfiles)

    def __exit__(self, exc_value: object, exc_type: object, traceback: object) -> None:
        self.stack.close()

def search(args):
    all_hmms = list(Path(args.markerset).iterdir())
    orthologs = {}
    seq_dict = {}
    for fasta in args.input:
        fasta = Path(fasta)
        with pyhmmer.easel.SequenceFile(fasta, digital=True) as seq_file, HMMFiles(*all_hmms) as hmm_file:
            sequences = seq_file.read_block()
            for seq in sequences[:]:
                seq_dict[seq.name] = seq.sequence
            for hits in pyhmmer.hmmsearch(hmm_file, sequences, E=args.evalue):
                cog = hits.query_name.decode()
                for hit in hits:
                    if hit.included:
                        if cog in orthologs:
                            orthologs[cog].add(hit.name.decode())
                        else:
                            orthologs[cog] = set([hit.name.decode()])
                        break
    
    
    count = 0
    for hits in orthologs.values():
        print(hits)
        if len(hits) >= 3:
            count += 1
            print(seq)
    print(f"Total orthologs found: {len(orthologs)}")
    print(f"Shared orthologs: {count}")
    print(f"seq_database: {len(seq_dict)}")
