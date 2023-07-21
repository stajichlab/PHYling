import logging
from pathlib import Path

import matplotlib.pyplot as plt

from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def tree_generator(file: Path, method: str) -> Phylo.BaseTree.Tree:
    MSA = AlignIO.read(file, format="fasta")
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, method)
    return constructor.build_tree(MSA)


def phylotree(inputs, input_dir, output, method, figure, **kwargs):
    """
    The tree module uses Biopython Phylo module to build a phylogenetic tree upon
    the multiple sequence alignemt (MSA) results.

    If passing multiple MSA results, the consensus tree method will be applied which
    use 50% cutoff to represent the majority of all the trees.

    By default will use UPGMA algorithm for tree building, you can switch to
    Neighbor Joining method by -m/--method nj as your preference.

    After building the tree, a ascii figure will be displayed and a treefile in
    newick format will be output. You can also use -f/--figure to output a
    matplotlib-style figure.
    """
    # If args.input_dir is used to instead of args.inputs
    logging.info(f"Algorithm choose for tree building: {method}")
    if input_dir:
        inputs = list(Path(input_dir).iterdir())
    else:
        inputs = [Path(sample) for sample in inputs]
    logging.info(f"Found {len(inputs)} MSA fasta")
    if inputs[0].name == "concat_alignments.faa":
        logging.info("Generate phylogenetic tree the on concatenated fasta")
    else:
        logging.info("Generate phylogenetic tree on all MSA fasta and conclude an majority consensus tree")
    output = Path(output)
    output.mkdir(exist_ok=True)

    if len(inputs) == 1:
        final_tree = tree_generator(inputs[0], method)
    else:
        trees = []
        for marker in inputs:
            trees.append(tree_generator(marker, method))
        final_tree = majority_consensus(trees, 0.5)

    Phylo.draw_ascii(final_tree)

    output_tree = output / f"{method}_tree.newick"
    logging.info(f"Output tree to {output_tree}")
    with open(output_tree, "w") as f:
        Phylo.write(final_tree, f, "newick")

    if figure:
        fig, ax = plt.subplots(figsize=(20, 12))
        output_fig = output / f"{method}_tree.png"
        Phylo.draw(final_tree, axes=ax)
        fig.savefig(output_fig)
