[![Conda build](https://img.shields.io/github/actions/workflow/status/stajichlab/PHYling/conda-build.yml?logo=github&label=conda%20build)](https://github.com/stajichlab/PHYling/actions/workflows/conda-build.yml)
[![Python build](https://img.shields.io/github/actions/workflow/status/stajichlab/PHYling/python-versions.yml?logo=github&label=python%20build)](https://github.com/stajichlab/PHYling/actions/workflows/python-versions.yml)
[![Python](https://img.shields.io/badge/python-3.9_%7C_3.10_%7C_3.11_%7C_3.12-blue?logo=python)](https://github.com/stajichlab/PHYling/actions/workflows/python-versions.yml)
[![codecov](https://codecov.io/gh/stajichlab/PHYling/graph/badge.svg?token=ZH5GBQYKZ6)](https://codecov.io/gh/stajichlab/PHYling)
[![License](https://img.shields.io/github/license/stajichlab/PHYling?label=license)](https://github.com/stajichlab/PHYling/blob/main/LICENSE)

# PHYling tool

The unified PHYling pipeline for phylogenomic data collection from annotated genomes.

This is latest iteration of tool for using phylogenetically conserved markers to pull out informative gene or protein info from
genomic and transcriptomic datasets in order to construct gene trees and species phylogenies.

The aligned markers can be extracted from protein sequences for phylogenetic analyses and also projected into coding sequence
alignments for codon-based analyses for better resolution of recently diverged species.

The assumptions in this approach are that the markers are generally single copy in genomes and taking best hit is sufficient first
approximation for identifying orthologs. A separate file is parsed and file best_multihits which lists all the hits above the
cutoff threshold for a given marker which can be used to assess duplication or attempt to incorporate paralogs into the analysis
down the road.

The marker sets developed for this approach in fungi are available as part of the [1KFG
Phylogenomics_HMMs](https://github.com/1KFG/Phylogenomics_HMMs) project resource and preferred use of the [BUSCO marker
sets][Busco].

### Flow chart

<img src= "misc/phyling_flowchart.svg" alt="PHYling flowchart" width="800">

### New features compared to the original version

- Using [pyhmmer] to improve the multithread performance in hmmsearch and hmmalign.
- Implement all stuff in python. The entire program will be more readable and maintainable.
- Simplify some steps and reduce the intermediate files as much as possible.
- [Muscle] is now available for alternative alignment method.
- Use [PhyKIT] to remove uninformative orthologs.
- [FastTree], [RAxML] and [IQTree] are now available for tree construction.
- [ASTER], a C++ version of [ASTRAL] is now integrated to resolve consensus among trees built upon individual genes.

## Usage

First of all, install the package following the [instruction](#install) below.

PHYling is a package to extract phylogenomic markers and build a phylogenetic tree upon them. It comprises 3 modules - download,
align and tree. Use `phyling --help` to see more details.

```
positional arguments:
  {download,align,tree}
    download            Download HMM markers
    align               Run multiple sequence alignments against orthologs found among samples
    tree                Build a phylogenetic tree based on multiple sequence alignment results

options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
```

To test run on the example files, please `cd` into the folder `example`.

```
cd example
```

In general, PHYling takes fasta as input. The gzipped fasta is also valid.

The folder `example/pep` includes 5 example peptide fasta which can be used for test run.

In addition to the peptide sequences, PHYling can also takes DNA coding sequences as inputs to more accurately estimate the
phylogeny of closely related species. When taking DNA coding sequences as inputs, DNA sequences will be translated into peptide
sequences and all the hmmsearch/align are done on the peptide version. The final MSA results will be back-translated into DNA at
the final stage. The example DNA sequences are placed under the folder `example/cds`.

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
```

Firstly, use `download list` to show the available BUSCO markersets.

```
phyling download list
```

By default the downloaded markersets will be saved to the `~/.phyling/HMM`. The **Datasets available online** section lists the
markersets that available on the [BUSCO][Busco] website.

And the **Datasets available on local** section lists the markersets that have already been downloaded.

To download the markerset, copy the name from the list and paste it to the download module directly. Here we use `fungi_odb10` as
example.

```
phyling download fungi_odb10
```

The download module will automatically check for updates to the markerset online each time it runs. The local markersets which
have available updates online will be marked as _\[Outdated\]_. You can rerun `phyling download [markerset]` to update the local
files.

### Find the orthologs and align them

The align module identify the orthologs among all the samples using _hmmsearch_. HMM profiles that have matches on more than 4
samples are considered **orthologs**.

Before conducting _hmmsearch_, the module will first search for the bitscore cutoff file within the root HMM folder. If the cutoff
file is not found, the reporting threshold for _hmmsearch_ will be determined based on the `-E/-evalue` (default is 1e-10).

Once the orthologs are identified, the sequences extracted from each sample undergo multiple sequence alignment. By default, the
alignment is performed using the _hmmalign_ method. However, users have the option to switch to _muscle_ by specifying the
`-M/--method muscle` flag.

Finally, each alignment result is output separately. You can decide whether you want to concatenate them or use consensus strategy
in the tree module. See all the options with `phyling align --help`.

```
options:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose mode for debug
  -i file [files ...], --inputs file [files ...]
                        Query pepetide/cds fasta or gzipped fasta
  -I directory, --input_dir directory
                        Directory containing query pepetide/cds fasta or gzipped fasta
  -o directory, --output directory
                        Output directory of the alignment results (default: phyling-align-20240423-173849-0700 [current timestamp])
  -m directory, --markerset directory
                        Directory of the HMM markerset
  -E float, --evalue float
                        Hmmsearch reporting threshold (default: 1e-10, only being used when bitscore cutoff file is not available)
  -M {hmmalign,muscle}, --method {hmmalign,muscle}
                        Program used for multiple sequence alignment (default: hmmalign)
  --non_trim            Report non-trimmed alignment results
  -t THREADS, --threads THREADS
                        Threads for hmmsearch and the number of parallelized jobs in MSA step. Better be multiple of 4 if using more than 8 threads (default: 1)
```

Run the align module with all the fasta files under folder `pep`.

```
phyling align -I pep -o align -m fungi_odb10
```

An equivalent way to send inputs.

```
phyling align -i pep/*.fasta.gz -o align -m fungi_odb10
```

Or if you're just interested in part of the fasta, you can specify the inputs one-by-one.

```
phyling align -i pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  pep/Zygorhynchus_heterogamous_NRRL_1489.aa.fasta.gz \
  -o align \
  -m fungi_odb10
```

**Note: Required at least 4 samples to build a tree!**

Accelerate by using 16 cpus.

```
phyling align -I pep -o align -m HMM/fungi_odb10/hmms -t 16
```

#### Multithreading strategy

According to [pyhmmer benchmark](https://pyhmmer.readthedocs.io/en/stable/benchmarks.html), the acceleration benefits from
multithreading drop significantly as more CPUs are utilized. When less then 8 cpus are given, the hmmsearch step will run on
single-thread manner and all cpus will be used for each round of hmmsearch. When 8 or more cpus are given, the hmmsearch step will
use 4 cpus for each parallel job. In the example above, 4 hmmsearch jobs will run parallelly and each job utilize 4 cpus. For the
alignment step, 16 parallel jobs will be launched and each parallel job is running on single-thread manner.

Highly recommended if **muscle** is chosen for alignment. (**muscle** is much slower than **hmmalign**!!)

#### Use coding sequences instead of peptide sequences

In some circumstances, the highly shared peptide sequences make it difficult to resolve the relationship among closely related
species. To address the issue, one can use DNA coding sequences (CDS), which contain more evolutionary traces, instead of peptide
sequences for phylogeny analysis.

Run the align module with cds fasta files under folder `cds`.

```
phyling align -I cds -o align_cds -m HMM/fungi_odb10/hmms -t 16
```

The CDS inputs will be translated into peptide sequences in the first steps. The translated peptide sequences will be used for
hmmsearch and the alignment steps. The peptide alignment results will then being back-translated according to the original CDS
inputs. And the back-translated DNA version alignments will be output.

#### Checkpoint for quick rerun

Once the align module complete, a checkpoint file will be generated to the output folder. This checkpoint file stores the
parameters, samples and identified orthologs, which will be loaded to the pipeline when rerun with the same output folder. Then
the align module will determine whether to skip the hmmsearch on some of the samples that were already completed in the previous
run.

For example, we run the align module by:

```
phyling align -i pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  pep/Zygorhynchus_heterogamous_NRRL_1489.aa.fasta.gz \
  -o align \
  -m fungi_odb10
  -t 16
```

And later we want to add `Actinomucor elegans` to the analysis but kick out `Zygorhynchus heterogamous`. We can start another run
and specifying the same output folder:

```
phyling align -i pep/Actinomucor_elegans_CBS_100.09.aa.fasta.gz
  pep/Pilobolus_umbonatus_NRRL_6349.aa.fasta.gz \
  pep/Rhizopus_homothallicus_CBS_336.62.aa.fasta.gz \
  pep/Rhizopus_rhizopodiformis_NRRL_2570.aa.fasta.gz \
  -o align \
  -m fungi_odb10
  -t 16
```

In this case, `Pilobolus umbonatus`, `Rhizopus homothallicus` and `Rhizopus rhizopodiformis` will be skipped from the hmmsearch
process since they have already been searched in the previous run. The `Actinomucor elegans` is the only sample need to be
hmmsearched. On the other hand, the `Zygorhynchus heterogamous` will be removed from the current run.

Note that if the input files have changes, they will also being detected by align module and trigger the rerun. If the hmmsearch
evalue/bitscore cutoff is changed, all the samples will need to rerun the hmmsearch step. Also, the changes on HMM markersets or
input samples with different seqtype will terminate the align module. (since this case should be considered an entirely different
analysis)

### Build tree from multiple sequence alignment results

Finally, we can run the tree module, use the multiple sequence alignment results to build a phylogenetic tree. By default, it uses
the _consensus tree_ strategy (conclude the majority of trees which was built upon each single gene) But you can choose to use
_concatenated alignment_ strategy by specifying `-c/--concat`. Currently, 3 methods (FastTree, RAxML and IQTree) are available for
tree building. You can choose your own preferred method by specifying `-M/--method`. (default is FastTree) We also use the
[treeness/RCV] calculated by [PhyKit] and select the most informative markers for final tree building. You can adjust the number
of selected markers by specifying `-n/--top_n_toverr`. See all the options with `phyling tree --help`.

```
options:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose mode for debug
  -i file [files ...], --inputs file [files ...]
                        Multiple sequence alignment fasta of the markers
  -I directory, --input_dir directory
                        Directory containing multiple sequence alignment fasta of the markers
  -o directory, --output directory
                        Output directory of the newick treefile (default: phyling-tree-20240423-183138-0700 [current timestamp])
  -M {ft,raxml,iqtree}, --method {ft,raxml,iqtree}
                        Algorithm used for tree building. (default: ft)
                        Available options:
                        FastTree: ft
                        RAxML: raxml
                        IQTree: iqtree
  -n TOP_N_TOVERR, --top_n_toverr TOP_N_TOVERR
                        Select the top n markers based on their treeness/RCV for final tree building (default: 50. Specify 0 to use all markers)
  -c, --concat          Concatenated alignment results
  -p {seq,codon,seq+codon}, --partition {seq,codon,seq+codon}
                        Create a partition file by sequence or by codon position when --concat enabled. "codon" and "seq+codon" only work when inputs are DNA sequences
  -f, --figure          Generate a matplotlib tree figure
  -t THREADS, --threads THREADS
                        Threads for tree construction (default: 1)
```

Run the tree module with all the alignment results under folder `align`.

```
phyling tree -I align
```

You can also use only part of the alignment results to build tree.

```
phyling tree -i align/100957at4751.aa.mfa align/174653at4751.aa.mfa align/255412at4751.aa.mfa
```

**Note: the tree module uses the checkpoint file .align.ckp in the input folder (or the parent folder of the input files) to
determine the sample names and seqtype. If the checkpoint file is missing or corrupted it can also automatically detects these
information but requires more time.**

Use IQTree instead of the default FastTree method for tree building and running with 16 threads.

```
phyling tree -I align -m iqtree -t 16
```

Use matplotlib to generate a tree figure.

```
phyling tree -I align -f -t 16
```

#### Use top_n_toverr to filter the markers by treeness/RCV

The align module sometimes reports a lot of orthologs depending on the size of BUSCO dataset and the gene set homogeneity among
samples. However, not every marker is equally informative in resolving phylogeny; some substitutions may not contribute
significantly to the branches in the phylogenetic tree. To address this issue, PHYling incorporates [PhyKIT] to compute the
[treeness/RCV] scores (toverr) and ranks the markers accordingly. Higher ranks indicate greater informativeness and lower
susceptibility to composition bias and will be selected for final tree building (thru _consensus_ or _concatenate_ strategy).

By default, PHYling pick the top 50 markers for further analysis but user can adjust the number by specifying `-n/--top_n_toverr`.
In the example below, we pick only the top 20 markers:

```
phyling tree -I align -n 20 -t 16
```

All markers will undergo tree building and compute their treeness/RCV scores, and the top 20 markers will be selected to
reconstruct the final consensus tree. (PHYling uses _consensus_ strategy by default)

Note that if the number of identified orthologs is less than the assigned `-n/--top_n_toverr` value, PHYling will use all markers
instead. Users can also specify 0 to use all the markers directly.

```
phyling tree -I align -n 0 -t 16
```

After tree construction, a file `top_toverr_trees.tsv` recording the treeness/RCV of selected markers and a folder `selected_MSAs`
containing the symlinks to the mfa of the selected markers will be output as well.

#### Use concatenate strategy

The consensus strategy is used for tree construction by default but users can choose to concatenate the markers and generate a
single tree on it.

```
phyling tree -I align -c -t 16
```

When concatenation mode is enabled, PHYling do the first round tree building on each marker with [FastTree] and compute their
[toverr][treeness/RCV]. These highly ranked markers are then selected for concatenation and do the final tree building. (using the
method specified by `-M/--method`) The concatnated fasta `concat_alignments.mfa` will also being output which allow users to
perform tree building with tools not incorporated in PHYling. The example below shows to pick the top 20 markers and build a tree
with _concatenate_ strategy by specifying `-n/--top_n_toverr`.

```
phyling tree -I align -n 20 -c -t 16
```

Meanwhile, users can construct tree with a more sophisticated **partition mode** when using [RAxML] and [IQTree]. In general, the
partition mode expects different sequence regions exhibit different evolutionary rates, which should be estimated with different
models. Here we provide 3 different most commonly-used modes:

- seq: partitioning by marker. (each marker evolves separately)
- codon: partitioning by codon of 3 of concatenated sequence. (only available on CDS. Each codon evolves separately)
- seq+codon: partitioning by codon of 3 marker-wisely. (only available on CDS. Each codon of each marker evolves separately)

The example below concatenate the top markers (50 by default) and run tree building with sequence partitioning thru iqtree.

```
phyling tree -I align -M iqtree -c -p seq -t 16
```

**Note: Partition mode is not supported in FastTree.**

#### Tune it yourselves

To adapt the most common needs, PHYling uses the very basic commands to run [FastTree], [RAxML] and [IQTree]:

FastTree for peptide:

```
FastTree -gamma file.mfa -lg
```

FastTree for cds:

```
FastTree -gamma file.mfa -nt
```

RAxML for peptide:

```
raxml -s file.mfa -p 12345 -w [absolute_path_output] -n pep -m PROTGAMMAAUTO
```

RAxML for peptide + partition mode:

```
raxml -s file.mfa -p 12345 -w [absolute_path_output] -n pep -m PROTGAMMAAUTO -q file.partition
```

RAxML for cds:

```
raxml -s file.mfa -p 12345 -w [absolute_path_output] -n cds -m GTRCAT
```

RAxML for cds + partition mode:

```
raxml -s file.mfa -p 12345 -w [absolute_path_output] -n cds -m GTRCAT -q file.partition
```

IQTree for pep/cds:

```
iqtree2 -s file.mfa --prefix [output_path] -m MFP --mset raxml -T AUTO
```

IQTree partition mode:

```
iqtree2 -s file.mfa --prefix [output_path] -m MFP --mset raxml -T AUTO -p file.partition
```

You can use the tree module to prepare the required data (i.e. concat_alignment.mfa or toverr filtered mfas) and rerun the tree
building step with your own preferred parameters.

## Requirements

- Python >= 3.9
- [Biopython](https://biopython.org/)
- [pyhmmer], a HMMER3 implementation on python3.
- [muscle] for alternative method for multiple sequence alignment. (Optional)
- [ClipKIT] for removing sites that are poor of phylogenetic signal.
- [PhyKIT] for calculating treeness/RCV to filter uninformative orthologs.
- [FastTree], use approximately maximum-likelihood to build trees.
- [RAxML], a more sophisticated maximum-likelihood-based tree building tool. (Optional)
- [IQTree], a modern maximum-likelihood-based tool for tree building. (Optional)
- [ASTER], a C++ re-implementation of [ASTRAL] to resolve consensus among trees.

## Install

Please download the source code from the [latest release](https://github.com/stajichlab/PHYling/releases/latest) and decompress it
or `git clone` the main branch.

Go into the PHYling folder. To avoid altering the base environment, it's advisable to install the software in a dedicated conda
environment Please use the environment.yml to create environment and install all the required packages.

```
cd PHYling-2.0.0
conda env create -f environment.yml
```

Install the package through pip in the PHYling folder.

```
pip install .
```

### Install additional package for developing (developer only)

Developer should clone the GitHub project directly instead of downloading from the releases. Some of the files for developing
purpose only are not included in the releases.

In addition to the requirements listed above, the following packages are required for developing environment.

- [pre-commit] >= 3.4.0

Developer can install it with conda by:

```
conda install pre-commit>=3.4.0
```

For convenience we also provide the additional conda env file. Please use the dev_additional_packages.yml to install the additional packages.

```
conda env update -f dev_additional_packages.yml
```

[Busco]: https://busco-data.ezlab.org/v5/data/lineages/
[pyhmmer]: https://pyhmmer.readthedocs.io/en/stable/index.html
[muscle]: https://drive5.com/muscle5/
[ClipKIT]: https://jlsteenwyk.com/ClipKIT/
[PhyKIT]: https://jlsteenwyk.com/PhyKIT/
[Treeness/RCV]: https://jlsteenwyk.com/PhyKIT/usage/index.html#treeness-over-rcv
[FastTree]: http://www.microbesonline.org/fasttree/
[RAxML]: https://cme.h-its.org/exelixis/web/software/raxml/
[IQTree]: http://www.iqtree.org/
[ASTER]: https://github.com/chaoszhang/ASTER
[ASTRAL]: https://github.com/smirarab/ASTRAL
[pre-commit]: https://pre-commit.com/
