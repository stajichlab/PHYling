"""Utilities that commonly used in different modules."""
import time
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


def runtime(start_time: float) -> str:
    """Measure the runtime of a process by given the start time."""
    elapsed_sec = time.monotonic() - start_time
    hours = int(elapsed_sec // 3600)
    mins = int((elapsed_sec // 60) % 60)
    secs = elapsed_sec % 60
    msg = f"{secs:.3f} seconds"
    msg = f"{mins} minutes " + msg if mins > 0 else msg
    msg = f"{hours} hours " + msg if hours > 0 else msg
    return msg
