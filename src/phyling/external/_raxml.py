"""RAxML utilities"""

import re
from pathlib import Path
from typing import Literal

from ..libphyling import TreeMethods
from ..libphyling._utils import check_binary
from ._base import TreeToolWrapper

RAXML_BIN = check_binary(
    TreeMethods.RAXML.method, TreeMethods.RAXML.bins, "bioconda::raxml-ng", "https://github.com/amkozlov/raxml-ng"
)


def _partition_raxml_to_nexus(file: str | Path, output: str | Path) -> str | Path:
    with open(file) as f:
        charset = []
        charpartition = []
        for line in f.read().strip().split("\n"):
            line = line.split(",")
            line[0] = re.sub(r"\+BU\{\d+\.\d+\}", "", line[0])
            line[0] = re.sub(r"\+FC|\+FO|\+FE|\+FU", "+F", line[0])
            line[0] = re.sub(r"\+IO|\+IC|\+IU", "+I", line[0])
            line[0] = re.sub(r"\+G4m", "+G4", line[0])
            line[0] = re.sub(r"\+GA", "+G", line[0])
            line[0] = re.sub(r"\{([^}]*)\}", lambda m: "{" + m.group(1).replace("/", ",") + "}", line[0])
            line = [i.strip() for i in [line[0], *line[1].split("=")]]
            charset.append(f"  charset {line[1]} = {line[2]}")
            charpartition.append(f"    {line[0]}: {line[1]}")
    if len(charpartition) == 1:
        return charpartition[0].split(":")[0].strip()
    with open(output, "w") as f:
        f.write("#nexus\nbegin sets;\n")
        f.write(";\n".join(charset) + ";\n")
        f.write("  charpartition mymodels =\n")
        f.write(",\n".join(charpartition) + ";\n")
        f.write("end;")
    return output


class Raxml(TreeToolWrapper):
    """Runs RAxML-NG to build a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree (MFA2Tree): The MFA2Tree object containing alignment data.
        output (str | Path | None, optional): Path to save the resulting tree file. If not provided, a temporary path is used.
        partition_file (str | Path | None, optional): Path to a partition file for model partitioning. Optional.
        bs (int): Bootstrap value. Defaults to 50.
        threads (int): Number of threads to use for RAxML-NG computation.
        capture_cmd (bool): If True, returns the RAxML-NG command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple of the Tree object and the command string.
    """

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        model: str = "AUTO",
        threads: int = 1,
        add_args: tuple | list | None = None,
    ):
        """Instantiate RAxML-NG runner.

        Args:
            file (str | Path): Path of the MSA file.
            output (str | Path): Directory to save the results.
            seqtype (str | Path | None, optional): The sequence type of the file.
            model (int): Model to use for phylogeny inference.
            threads (int): Number of threads to use.
            add_args (bool): Additional arguments to pass in.

        Returns:
            If capture_cmd is False (default), returns a Tree object.
            If capture_cmd is True, returns a tuple of the Tree object and the command string.
        """
        super().__init__(TreeMethods.RAXML.method, file, output, seqtype=seqtype, model=model, add_args=add_args, threads=threads)

    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        *,
        seqtype: Literal["DNA", "AA"] | None,
        model: str = "AUTO",
        threads: int = 1,
    ):
        self._cmd = [
            RAXML_BIN,
            "--msa",
            str(file),
            "--prefix",
            str(output),
            "--model",
            str(model),
            "--threads",
            f"auto{{{threads}}}",
        ]
        if seqtype:
            self._cmd.extend(["--data-type", seqtype])

        self._target = Path(f"{output}.raxml.bestTree")

    def _update_model(self):
        model_file = self._target.with_suffix(".bestModel")
        self._model = _partition_raxml_to_nexus(model_file, f"{model_file}.nex")
