"""PHYling - Phylogenomic reconstruction from genomes.

PHYling comprises 4 modules - download, align, filter and tree. The download module can be used to download HMM markerset from
BUSCO. The align module is the core element of this package which generate multiple sequence alignment among the orthologs
found across samples. The filter module calculates treeness/RCV scores to filter out the uninformative alignment results. The
tree module help to build a phylogenetic tree by different algorithms.
"""

import logging
import os
from importlib.metadata import metadata
from pathlib import Path

# Create logger for the package
logger = logging.basicConfig(format="%(asctime)s %(name)s %(levelname)s %(message)s", level="INFO")
logger = logging.getLogger(__name__)


VERSION = metadata("phyling")["Version"]
AUTHOR = metadata("phyling")["Author-email"]

# Create config folder in $HOME/.phyling
CFG_DIRS: list[Path] = []
if os.getenv("PHYLING_DB"):
    CFG_DIRS.extend([Path(path) for path in os.getenv("PHYLING_DB").split(":")])
if CFG_DIRS and os.access(CFG_DIRS[0], os.W_OK):
    pass
else:
    cfg_dir = Path.home() / ".phyling"
    cfg_dir.mkdir(exist_ok=True)
    CFG_DIRS.insert(0, cfg_dir)

# Available CPUs
AVAIL_CPUS = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count()))
