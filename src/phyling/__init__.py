"""Placeholder init file for Python purposes."""
from importlib.metadata import version

import phyling.config as config

__name__ = "PHYling"
__version__ = version("phyling")

# Create config folder in $HOME/.phyling
config.cfg_dir.mkdir(exist_ok=True)
