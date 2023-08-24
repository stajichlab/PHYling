![Conda](https://github.com/stajichlab/PHYling_unified/actions/workflows/conda-building-check.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.9_%7C_3.10_%7C_3.11-blue)

# PHYling tool
The unified PHYling pipeline for phylogenomic data collection from annotated genomes.

This is latest iteration of tool for using phylogenetically conserved markers to pull out informative
gene or protein info from genomic and transcriptomic datasets in order to construct gene trees and species phylogenies.

The aligned markers can be extracted from protein sequences for phylogenetic analyses and also projected into coding sequence alignments for codon-based analyses for better resolution of recently diverged species.

The assumptions in this approach are that the markers are generally single copy in genomes and taking best hit is sufficient first approximation for identifying orthologs. A separate file is parsed and file best_multihits which lists all the hits above the cutoff threshold for a given marker which can be used to assess duplication or attempt to incorporate paralogs into the analysis down the road.

The marker sets developed for this approach in fungi are available as part of the [1KFG Phylogenomics_HMMs](https://github.com/1KFG/Phylogenomics_HMMs) project resource and preferred use of the [BUSCO marker sets][Busco].

### Flow chart
![PHYling flowchart](misc/phyling_flowchart.png)

### New features compared to the original version
- Using [pyhmmer] to improve the multithread performance in hmmsearch and hmmalign.
- Implement all stuff in python. The entire program will be more readable and maintainable.
- Simplify some steps and reduce the intermediate files as much as possible.
- [Muscle] is now available for alternative alignment method.
- [VeryFastTree] is now available for tree construction.
- [ASTER], a C++ version of [ASTRAL] is now integrated to resolve consensus among trees built upon individual genes.

## Usage
First of all, install the package following the [instruction](#requirements-and-installation) below.

PHYling is a package to extract phylogenomic markers and build a phylogenetic
tree upon them. It comprises 3 modules - download, align and tree. Use `phyling --help` to see more details.
```
positional arguments:
  {download,align,tree}
    download            Download HMM markers
    align               Run multiple sequence alignments against orthologs found among samples
    tree                Build a phylogenetic tree based on multiple sequence alignment results

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
```

To test run on the example files, please `cd` into the folder `example`.
```
cd example
```
In general, PHYling takes fasta as input. The gzipped fasta is also valid.

All the search and alignment processes are done on peptide sequences.
The folder `example/pep` includes 5 example peptide fasta which can be used for test run.

In addition to the peptide sequences, PHYling can also takes DNA coding sequences as inputs to more accurately estimate the phylogeny of closely related species.
If you're interested in testing the feature, there are also some example DNA sequences under the folder `example/cds`.

### Download HMM markerset
The download module is used to download the HMM markerset from BUSCO website. (Currently is updated to [v5][Busco])
See all options with `phyling download --help`.
```
positional arguments:
  HMM markerset or "list"
                        Name of the HMM markerset

options:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose mode for debug
  -o OUTPUT, --output OUTPUT
                        Output directory to save HMM markerset (default="./HMM")
```

Firstly, use `download list` to show the available BUSCO markersets.
```bash
phyling download list
```

Copy and paste the markerset to download it. Here we use `fungi_odb10` as example.
```bash
phyling download fungi_odb10
```

### Find the orthologs and align them
The align module identify the orthologs among all the samples using *hmmsearch*.
HMM profiles that have matches on more than 3 samples are considered **orthologs**.

Before conducting *hmmsearch*, the module will first search for the bitscore cutoff file within the root HMM folder.
If the cutoff file is not found, the reporting threshold for *hmmsearch* will be determined based on the `-E/-evalue` (default is 1e-10).

Once the orthologs are identified, the sequences extracted from each sample undergo multiple sequence alignment.
By default, the alignment is performed using the *hmmalign* method.
However, users have the option to switch to using *muscle* by specifying the `-M/--method muscle` flag.

By default, each alignment result is output separately and is expected to resolve their phylogeny by consensus tree method.
If you prefer to use concatenate strategy. You can concatenate all the alignment by passing `-c/--concat`.
See all the options with `phyling align --help`.
```
options:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose mode for debug
  -i INPUTS [INPUTS ...], --inputs INPUTS [INPUTS ...]
                        Query pepetide/cds fasta or gzipped fasta
  -I INPUT_DIR, --input_dir INPUT_DIR
                        Directory containing query pepetide/cds fasta or gzipped fasta
  -o OUTPUT, --output OUTPUT
                        Output directory of the alignment results (default="./align")
  -m MARKERSET, --markerset MARKERSET
                        Directory of the HMM markerset
  -E EVALUE, --evalue EVALUE
                        Hmmsearch reporting threshold (default=1e-10)
  -M {hmmalign,muscle}, --method {hmmalign,muscle}
                        Program used for multiple sequence alignment (default="hmmalign")
  -n, --non_trim        Report non-clipkit-trimmed alignment results
  -c, --concat          Report concatenated alignment results
  -t THREADS, --threads THREADS
                        Threads for hmmsearch and the number of parallelized jobs in MSA step (default=1)
```

Run the align module with all the fasta files under folder `pep`.
```bash
phyling align -I pep -m HMM/fungi_odb10/hmms
```

An equivalent way to send inputs.
```bash
phyling align -i pep/*.fasta.gz -m HMM/fungi_odb10/hmms
```

Or if you're just interested in part of the fasta, you can specify the inputs one-by-one.
```bash
phyling align -i pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  -m HMM/fungi_odb10/hmms
```
**Note: Required at least 3 samples in order to build a tree!**

Accelerate by using 16 cpus.
```bash
phyling align -I pep -m HMM/fungi_odb10/hmms -t 16
```
According to [pyhmmer benchmark](https://pyhmmer.readthedocs.io/en/stable/benchmarks.html) the multithread acceleration drop dramatically when more cpu is used.
When less then 8 cpus are given, the hmmsearch step will run in single-thread manner and all cpus will be used for each round of hmmsearch.
When 8 or more cpus are given, the hmmsearch step will use 4 cpus for each parallel job.
In this example, 4 hmmsearch jobs will run parallelly and each job utilize 4 cpus.
For the alignment step, 16 parallel jobs will be launched and each parallel job is running in single-thread manner.

Highly recommended if **muscle** is chosen for alignment. (**muscle** is much slower than **hmmalign**!!)

#### Use coding sequence instead of peptide sequence
In some circumstances, the highly shared peptide sequences make it difficult to resolve the relationship among closely related species.
To address the issue, one can use DNA coding sequences (CDS), which contain more evolutionary traces, instead of peptide sequences for phylogeny analysis.

Run the align module with cds fasta files under folder `cds`.
```bash
phyling align -I cds -m HMM/fungi_odb10/hmms -t 16
```
The CDS inputs will be translated into peptide sequences in the first steps.
The translated peptide sequences will be used for hmmsearch and the alignment steps.
The peptide alignment results will then being back translated according to the original CDS inputs.
And the back-translated DNA version alignments will be output.

### Build tree on multiple sequence alignment results
Finally, we can run the tree module to use the multiple sequence alignment results to build a phylogenetic tree.
It support both *consensus tree* (conclude the majority of trees which was built upon each single gene) and *concatenated alignment* method.
See all the options with `phyling tree --help`.
```
options:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose mode for debug
  -i INPUTS [INPUTS ...], --inputs INPUTS [INPUTS ...]
                        Multiple sequence alignment fasta
  -I INPUT_DIR, --input_dir INPUT_DIR
                        Directory containing multiple sequence alignment fasta
  -o OUTPUT, --output OUTPUT
                        Output directory of the newick treefile (default=".")
  -M {upgma,nj,ft}, --method {upgma,nj,ft}
                        Algorithm used for tree building (default="upgma")
  -f, --figure          Generate a matplotlib tree figure
  -t THREADS, --threads THREADS
                        Threads for tree construction (default=1)
```

Run the tree module with all the alignment results under folder `align`.
```bash
phyling tree -I align
```

You can also use only part of the alignment results to build tree.
```bash
phyling tree -i align/100957at4751.faa align/174653at4751.faa align/255412at4751.faa
```

Use VeryFastTree instead of the default UPGMA method for tree building and running with 16 threads.
```bash
phyling tree -I align -m ft -t 16
```

Use matplotlib to generate a tree figure.
```bash
phyling tree -I align -f
```

## Requirements and Installation
- Python >= 3.9
- [Biopython](https://biopython.org/)
- [pyhmmer], a HMMER3 implementation on python3.
- [clipkit](https://jlsteenwyk.com/ClipKIT/) for trimming.
- [muscle] for alternative method for multiple sequence alignment. (Optional)
- [VeryFastTree], use approximately maximum-likelihood to build trees. (Optional)
- [ASTER], a C++ re-implementation of ASTRAL(https://github.com/smirarab/ASTRAL) to resolve consensus among trees.

Use the environment.yml to install all the required packages
```bash
conda env create -f environment.yml
```

Go into the phyling folder and install the package through pip
```bash
pip install .
```

## Notes
- Training your own marker set is also possible but most busco sets are good starting place.
- The multiple sequence alignment results can also be sent to other phylogenetic tool like IQ-tree for tree building.

[Busco]: https://busco-data.ezlab.org/v5/data/lineages/
[pyhmmer]: https://pyhmmer.readthedocs.io/en/stable/index.html
[muscle]: https://drive5.com/muscle5/
[VeryFastTree]: https://github.com/citiususc/veryfasttree
[ASTER]: https://github.com/chaoszhang/ASTER
