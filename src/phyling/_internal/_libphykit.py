"""Phykit utilities."""

from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import itertools
from pathlib import Path

import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from phykit.services.alignment.base import Alignment as Phykit_alignment
from phykit.services.tree import Saturation as Phykit_saturation
from phykit.services.tree.base import Tree as Phykit_tree
from sklearn.linear_model import LinearRegression


class Saturation(Phykit_saturation):
    """Class for calculating saturation."""

    def __init__(self) -> None:
        """Initiate the object."""
        pass

    def compute_saturation(self, alignment: MultipleSeqAlignment, tree: Tree) -> float:
        """Calculate the saturation of the tree by PhyKIT implementation."""
        tips = Phykit_tree().get_tip_names_from_tree(tree=tree)
        combos = list(itertools.combinations(tips, 2))

        patristic_distances, uncorrected_distances = self.loop_through_combos_and_calculate_pds_and_pis(combos, alignment, tree)

        model = LinearRegression(fit_intercept=False)
        model.fit(np.array(patristic_distances).reshape(-1, 1), np.array(uncorrected_distances))
        slope = model.coef_[0]

        return round(slope, 4)


def compute_toverr(alignment: Path, tree: Tree) -> float:
    """Calculate the treeness/RCV of the tree by PhyKIT implementation."""
    # calculate treeness
    treeness = Phykit_tree().calculate_treeness(tree=tree)

    # calculate rcv
    aln = Phykit_alignment(alignment_file_path=alignment)
    relative_composition_variability = aln.calculate_rcv()

    # calculate treeness/rcv
    return round(treeness / relative_composition_variability, 4)
