"""Supporting routines for downloading data especially BUSCO markers."""

from __future__ import annotations

import csv
import hashlib
import shutil
import tarfile
import tempfile
from contextlib import ContextDecorator
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import phyling._internal._config as _config
from phyling import logger


class BuscoMetadata(ContextDecorator):
    """Store/update local copy of HMM files into destination."""

    def __init__(self, database_url: str, metadata: str | Path) -> None:
        """Initiate the object and download the latest metadata from online database."""
        self._hmm_folder = _config.cfg_dir / _config.default_HMM
        self._metadata_file = Path(metadata)
        self._get_metadata_online(database_url)

    def __enter__(self) -> BuscoMetadata:
        """Define the actions that will run when the object is created with `with` statement."""
        self._get_local_metadata()
        return self

    def __exit__(self, *exc) -> None:
        """Define the actions that will run when the object is going to be destroyed by the end of the `with` statement."""
        with open(self._metadata_file, "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["#name", "path", "md5"])
            for name, info in self._local_metadata.items():
                writer.writerow([name, info["path"], info["md5"]])

    @property
    def online(self) -> list[str]:
        """Return a list of the markersets retrieve from online metadata."""
        return sorted(list(self._online_metadata))

    @property
    def local(self) -> list[str]:
        """Return a list of the markersets retrieve from local metadata."""
        local_markersets = []
        for markerset, info in self._local_metadata.items():
            if info["md5"] != self._online_metadata[markerset]["md5"]:
                markerset += " [Outdated]"
            local_markersets.append(markerset)
        return sorted(local_markersets)

    def download(self, markerset: str) -> None:
        """Download the given markerset from busco database."""
        if markerset in self._local_metadata:
            if self._local_metadata[markerset]["md5"] == self._online_metadata[markerset]["md5"]:
                raise FileExistsError("Markerset already exists and update to date.")
            logger.info("Local markerset outdated. Update from online database...")
        logger.debug("Markerset not found. Download from online database...")
        data = fetch_url(self._online_metadata[markerset]["url"])
        md5 = hashlib.md5(data).hexdigest()
        if md5 != self._online_metadata[markerset]["md5"]:
            raise HTTPError("The checksum of the downloaded content doesn't match to the metadata.")
        file_path = self._save_markerset(data, markerset=markerset)
        self._local_metadata[markerset] = {"path": file_path, "md5": md5}

    def _save_markerset(self, data: bytes, markerset: str) -> Path:
        """Save the content of the markerset to a local folder."""
        output = self._hmm_folder / markerset
        if output.exists() and output.is_dir():
            shutil.rmtree(output)
        output.mkdir(parents=True)
        with tempfile.NamedTemporaryFile("wb", delete=False) as fp:
            fp.write(data)
            fp.close()
            logger.debug(f"Extract files to {output}.")
            with tarfile.open(fp.name, "r:gz") as f:
                f.extractall(output.parent, filter="data")
        Path(fp.name).unlink()
        return output

    def _get_metadata_online(self, busco_database_url: str):
        """Get the metadata from busco url."""
        self._online_metadata = {}
        data = fetch_url(f"{busco_database_url}/file_versions.tsv")
        for line in data.decode().split("\n"):
            line = line.split("\t")
            if line[-1] == "lineages":
                self._online_metadata[line[0]] = {
                    "url": f"{_config.database}/lineages/{line[0]}.{line[1]}.tar.gz",
                    "md5": line[2],
                }

    def _get_local_metadata(self):
        """Get the metadata from local file."""
        self._local_metadata = {}
        if self._metadata_file.exists() and self._metadata_file.is_file():
            with open(self._metadata_file) as f:
                for line in csv.reader(f, delimiter="\t"):
                    if line[0].startswith("#"):
                        continue
                    if not (self._hmm_folder / line[0]).exists():
                        continue
                    self._local_metadata[line[0]] = {"path": line[1], "md5": line[2]}
        else:
            logger.debug("Local metadata not found. Please download from the online database.")


def fetch_url(url: str) -> bytes:
    """Fetch URL data content.

    Attributes
    ----------
    url : str
        The url source.

    Return
    ------
    bytes
    """
    try:
        logger.debug(f"Download from {url} ...")
        with urlopen(url) as response:
            content = response.read()
    except HTTPError as e:
        raise HTTPError("URL currently unavailble.") from e
    except URLError as e:
        raise URLError("URL not found or no internet connection available.") from e
    return content


def wrapper(item_list: list[str], col: int, col_width: int, msg: str | None = None) -> None:
    """Adjust databases display according to the terminal size."""
    item_list = [item_list[x : x + col] for x in range(0, len(item_list), col)]
    if msg:
        print(msg)
        print()
    for row in item_list:
        # Print the database list
        print(" ".join(word.ljust(col_width) for word in row))
    print()
