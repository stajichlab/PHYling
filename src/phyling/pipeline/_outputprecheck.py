"""Output precheck subclasses for align and tree modules."""

from __future__ import annotations

from pathlib import Path

from .. import logger
from ..libphyling import TreeOutputFiles, _abc
from ..libphyling._utils import remove_dirs, remove_files
from ..libphyling.align import SampleList, SampleSeqs, SearchHitsManager
from ..libphyling.tree import MFA2Tree, MFA2TreeList


class AlignPrecheck(_abc.OutputPrecheckABC):
    """A single-instance class that provides features for input/output precheck, checkpoint loading/saving and final MSA output.

    Attributes:
        output (Path): The path to the output directory.
        ckp (str): The file name of the checkpoint.
        params (dict): A dictionary of parameters.
        samplelist (SampleList): List of sample sequences for processing.
    """

    __slots__ = ("samplelist",)

    ckp: str = ".align.ckp"
    params: dict = {param: None for param in ("markerset", "markerset_cutoff", "method")}

    def __init__(self, output: Path, samplelist: SampleList, **params) -> None:
        """Initialize AlignPrecheck with output directory and SampleList.

        Args:
            output (Path): Directory for storing output files.
            samplelist (SampleList): A collection of sample sequences.
            **params: Additional parameters for alignment.

        Raises:
            TypeError: If `samplelist` is not an instance of `SampleList`.
        """
        super().__init__(output)
        if not isinstance(samplelist, SampleList):
            raise TypeError(f"Argument samplelist only accepts {SampleList.__qualname__}. Got {type(samplelist)}.")
        self.samplelist: SampleList[SampleSeqs] = samplelist
        self.params.update(**params)

    def precheck(self, *, force_rerun: bool = False) -> tuple[SampleList, SearchHitsManager]:
        """Check the output folder and decide rerun actions for orthologs and samples.

        Args:
            force_rerun (bool, optional): Force rerun regardless of checkpoint status. Defaults to False.

        Returns:
            tuple[SampleList, SearchHitsManager]: Updated sample list and search hits.
        """
        super().precheck(force_rerun=force_rerun)
        if not any(self.output.iterdir()):
            return self.samplelist, SearchHitsManager()
        prev_params, prev_searchhits = self.load_checkpoint()
        return self._determine_rerun(prev_params, prev_searchhits)

    def load_checkpoint(self) -> tuple[dict, SearchHitsManager]:
        """Load checkpoint data for rerun decisions.

        Returns:
            tuple[dict, SearchHitsManager]: Loaded parameters and search hits.

        Raises:
            RuntimeError: If the checkpoint file is corrupted.
        """
        try:
            prev_params, prev_searchhits = super().load_checkpoint()
        except RuntimeError:
            raise RuntimeError(
                f"Checkpoint file {self.output / self.ckp} is corrupted. Please remove the output folder and try again."
            )
        self._type_check(prev_params, prev_searchhits)
        return prev_params, prev_searchhits

    def save_checkpoint(self, searchhits: SearchHitsManager) -> None:
        """Save the current parameters and search hits as a checkpoint.

        Args:
            searchhits (SearchHitsManager): Search hits data to save.

        Raises:
            RuntimeError: If the checkpoint data fails type validation.
        """
        self._type_check(self.params, searchhits)
        return super().save_checkpoint(searchhits)

    def _type_check(self, params: dict, searchhits: SearchHitsManager) -> None:
        """Validate the types of parameters and search hits.

        Args:
            params (dict): Parameters to validate.
            searchhits (SearchHitsManager): Search hits to validate.

        Raises:
            RuntimeError: If the types are incorrect.
        """
        if not isinstance(params, dict) or not isinstance(searchhits, SearchHitsManager):
            raise RuntimeError(
                f"Checkpoint file {self.output / self.ckp} is corrupted. Please remove the output folder and try again."
            )

    def _determine_rerun(
        self,
        prev_params: dict,
        prev_searchhits: SearchHitsManager,
    ) -> tuple[SampleList, SearchHitsManager]:
        """Handle rerun logic based on checkpoint data.

        Args:
            prev_params (dict): Previously saved parameters.
            prev_searchhits (SearchHitsManager): Previously saved search hits.

        Returns:
            tuple[SampleList, SearchHitsManager]: Updated sample list and filtered search hits.

        Raises:
            SystemExit: If incompatible changes are detected in the parameters or sequence type.
        """
        prev_samplelist = prev_searchhits.samplelist

        if prev_samplelist.seqtype != self.samplelist.seqtype:
            raise SystemExit("Seqtype is changed. Aborted.")
        diff_params = {param[0] for param in set(self.params.items()) ^ set(prev_params.items())}
        diff_inputs = set(prev_samplelist.checksums.items()) ^ set(self.samplelist.checksums.items())
        if not (diff_params or diff_inputs) and sum(1 for _ in self.output.iterdir()) > 1:
            raise SystemExit("Files not changed and parameters are identical to the previous run. Aborted.")
        if "markerset" in diff_params:
            raise SystemExit("Markerset is changed. Aborted.")
        if "markerset_cutoff" in diff_params:
            return self.samplelist, SearchHitsManager()

        # Remove files
        rm_files = [x for x in self.output.iterdir() if x.is_file() if x.name != self.ckp]
        drop_samples = [x.name for x in prev_samplelist if x.name not in self.samplelist]
        if drop_samples:
            logger.info("Remove samples that no longer exist from the orthologs collection retrieved from the checkpoint.")
            logger.debug("Removed samples: %s", ", ".join(drop_samples))
            searchhits = prev_searchhits.filter(drop_samples=drop_samples)
        else:
            searchhits = prev_searchhits
        remove_files(*rm_files)

        # Determine the samples need rerun
        remained_samples = SampleList()
        for sample in self.samplelist:
            if sample in searchhits.samplelist:
                if sample.checksum == searchhits.samplelist[sample.name].checksum:
                    continue
            remained_samples.append(sample)

        return remained_samples, searchhits


class FilterPrecheck(_abc.OutputPrecheckABC):
    """A class providing input/output precheck, checkpoint management, and filtered MSAs output.

    Attributes:
        output (Path): Path to the output directory.
        ckp (str): Filename for checkpoint storage.
        params (dict): Dictionary of parameters related to tree building.
        mfa2treelist (MFA2TreeList): List of MFA2Tree objects for processing.
    """

    ckp: str = ".filter.ckp"
    params: dict = {"top_n_toverr": None}

    def __init__(self, output: Path, mfa2treelist: MFA2TreeList, **params) -> None:
        """Initialize TreePrecheck with output directory and MFA2Tree list.

        Args:
            output (Path): Directory for storing output files.
            mfa2treelist (MFA2TreeList): A collection of MFA2Tree objects.
            **params: Additional parameters for tree generation.

        Raises:
            TypeError: If `mfa2treelist` is not an instance of `MFA2TreeList`.
        """
        super().__init__(output)
        if not isinstance(mfa2treelist, MFA2TreeList):
            raise TypeError(f"Argument mfa2treelist only accepts {MFA2TreeList.__qualname__}. Got {type(mfa2treelist)}.")
        self.mfa2treelist: MFA2TreeList[MFA2Tree] = mfa2treelist
        self.params.update(**params)

    def precheck(self, *, force_rerun: bool = False) -> tuple[MFA2TreeList, MFA2TreeList]:
        """Check the output folder and determine files and trees to rerun.

        Args:
            force_rerun (bool, optional): Force rerun regardless of checkpoint status. Defaults to False.

        Returns:
            tuple[MFA2TreeList, MFA2TreeList]: Updated MFA2Tree list and completed trees.
        """
        super().precheck(force_rerun=force_rerun)
        if not any(self.output.iterdir()):
            return self.mfa2treelist, MFA2TreeList()
        prev_params, prev_mfa2treelist = self.load_checkpoint()
        return self._determine_rerun(prev_params, prev_mfa2treelist)

    def load_checkpoint(self) -> tuple[dict, MFA2TreeList]:
        """Load checkpoint data for rerun decisions.

        Returns:
            tuple[dict, MFA2TreeList]: Loaded parameters and MFA2Tree list.

        Raises:
            RuntimeError: If the checkpoint file is corrupted.
        """
        try:
            prev_params, mfa2treelist = super().load_checkpoint()
        except RuntimeError:
            raise RuntimeError(
                f"Checkpoint file {self.output / self.ckp} is corrupted. Please remove the output folder and try again."
            )
        self._type_check(prev_params, mfa2treelist)
        return prev_params, mfa2treelist

    def save_checkpoint(self, mfa2treelist: MFA2TreeList) -> None:
        """Save the current parameters and MFA2Tree list as a checkpoint.

        Args:
            mfa2treelist (MFA2TreeList): MFA2Tree list data to save.

        Raises:
            RuntimeError: If the checkpoint data fails type validation.
        """
        self._type_check(self.params, mfa2treelist)
        super().save_checkpoint(mfa2treelist)

    def _type_check(self, params, mfa2treelist):
        """Validate the types of parameters and MFA2Tree list.

        Args:
            params (dict): Parameters to validate.
            mfa2treelist (MFA2TreeList): MFA2Tree list to validate.

        Raises:
            RuntimeError: If the types are incorrect.
        """
        if not isinstance(params, dict) or not isinstance(mfa2treelist, MFA2TreeList):
            raise RuntimeError(
                f"Checkpoint file {self.output / self.ckp} is corrupted. Please remove the output folder and try again."
            )

    def _determine_rerun(self, prev_params: dict, prev_mfa2treelist: MFA2TreeList) -> MFA2TreeList:
        """Handle rerun logic based on checkpoint data.

        Args:
            prev_params (dict): Previously saved parameters.
            prev_mfa2treelist (MFA2TreeList): Previously saved MFA2Tree list.

        Returns:
            tuple[MFA2TreeList, MFA2TreeList]: Updated MFA2Tree list and completed trees.

        Raises:
            SystemExit: If incompatible changes are detected in the parameters or sequence type.
        """

        if prev_mfa2treelist.seqtype != self.mfa2treelist.seqtype:
            raise SystemExit("Seqtype is changed. Aborted.")
        diff_params = {param[0] for param in set(self.params.items()) ^ set(prev_params.items())}
        if not diff_params and len(list(self.output.iterdir())) > 1:
            raise SystemExit("Files not changed and parameters are identical to the previous run. Aborted.")

        completed_mfa2treelist = MFA2TreeList()
        remained_mfa2treelist = MFA2TreeList()
        for msa in self.mfa2treelist:
            if msa.name in prev_mfa2treelist:
                prev_msa = prev_mfa2treelist[msa.name]
                if msa.checksum == prev_msa.checksum and prev_msa.toverr:
                    completed_mfa2treelist.append(prev_msa)
                    continue
            remained_mfa2treelist.append(msa)

        # Remove files
        rm_files = [x for x in self.output.iterdir() if x.is_file() if x.name == TreeOutputFiles.TREENESS]
        rm_dirs = [x for x in self.output.iterdir() if x.is_dir if x.name == TreeOutputFiles.MSAS_DIR]
        remove_files(*rm_files)
        remove_dirs(*rm_dirs)

        return remained_mfa2treelist, completed_mfa2treelist
