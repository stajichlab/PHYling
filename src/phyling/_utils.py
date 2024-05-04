"""Utilities that commonly used in different modules."""
from __future__ import annotations

import logging
import time
from functools import wraps
from hashlib import md5
from pathlib import Path
from typing import Iterable
from zlib import crc32


def get_file_checksum(file: str | Path) -> int:
    """Compute the crc checksum of a file."""
    file = Path(file)
    if not file.exists():
        raise FileNotFoundError(f"{file} not found.")
    crc = 0
    with open(file, "rb") as f:
        crc += crc32(f.read())
    return crc


def get_multifiles_checksum(files: Iterable[str | Path]) -> str:
    """
    Compute the md5 checksum of the given files.

    Crc checksum will be generated for each file and being sum up to compute the representative md5 checksum.
    """
    crcsum = 0
    for file in files:
        crcsum += get_file_checksum(file)
    return md5(str(crcsum).encode()).hexdigest()


def timing(func):
    """A decorator to measure the elapsed time of a function."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        """The wrapper function that measure the elapsed time."""
        t_start = time.monotonic()
        result = func(*args, **kwargs)
        t_end = time.monotonic()
        elapsed_sec = t_end - t_start
        hours = int(elapsed_sec // 3600)
        mins = int((elapsed_sec // 60) % 60)
        secs = elapsed_sec % 60
        msg = f"{secs:.3f} seconds"
        msg = f"{mins} minutes " + msg if mins > 0 else msg
        msg = f"{hours} hours " + msg if hours > 0 else msg
        logging.debug(f"{func.__module__}.{func.__name__} finished in {msg}.")
        return result

    return wrapper


def check_cls_vars(*vars: str):
    """A decorator to check whether a class variable is already defined before calling a class method."""

    def decorator(func):
        """The actual decorator function that wraps the class method."""

        @wraps(func)
        def wrapper(cls: type, *args, **kwargs):
            """The wrapper function that checks the class variable and calls the method."""
            for var in vars:
                if not hasattr(cls, var):
                    raise ValueError(
                        f'The class variable "{var}" of "{cls.__module__}.{cls.__name__}" must be set before calling this method.'
                    )
            return func(cls, *args, **kwargs)

        return wrapper

    return decorator
