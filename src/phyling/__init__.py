"""Required actions when initiating the package."""

import logging
from importlib.metadata import metadata

import phyling._internal._config as _config

_meta = metadata("phyling")
__version__ = _meta["Version"]
__author__ = _meta["Author-email"]
__all__ = []

# Create logger for the package
logger = logging.basicConfig(format="%(asctime)s %(name)s %(levelname)s %(message)s", level="INFO")
logger = logging.getLogger(__name__)

# Create config folder in $HOME/.phyling
_config.cfg_dir.mkdir(exist_ok=True)
