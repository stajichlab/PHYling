"""Supporting routines for downloading data especially BUSCO markers."""
from __future__ import annotations

import hashlib
import logging
import pickle
import shutil
import sys
import tarfile
from abc import ABC, abstractmethod
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

import phyling.config


class Data_updater(ABC):
    """Store, update and retrieve BUSCO markers."""

    def __init__(self, database_url, **kwargs):
        """Initialize database URL."""
        self._database_url = database_url

    @property
    @abstractmethod
    def _get_local_md5(self):
        """Compute md5 of local downloads to check as to whether to update."""
        pass

    @property
    @abstractmethod
    def _get_remote_md5(self):
        """Check remote md5 for downloads to check as to whether to update."""
        pass

    def _load_data(self):
        """Load datasets up."""
        pass

    @abstractmethod
    def _save_data(self, data):
        """Save datasets when downloading."""
        pass

    def fetch_url(self, url) -> bytes:
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
            logging.debug(f"Download from {url} ...")
            with urlopen(url) as response:
                content = response.read()
        except HTTPError:
            logging.error("URL not found or currently unavailble")
            sys.exit(1)
        return content

    def updater(self) -> dict:
        """Determine whether to update local copy of data download."""
        logging.debug(f"Check {self._filetype} exist or not ...")
        if self._data.is_file():
            logging.debug(f"{self._filetype} exist ({self._data})")
            logging.debug(f"Compare md5 checksums between the online latest record and the local {self._filetype}")
            if self._get_local_md5 == self._get_remote_md5:
                logging.debug("md5 is the same. No need to update metadata")
                return self._load_data()
            else:
                logging.info(
                    f"The local {self._filetype} md5 is different from the online latest record. "
                    f"Update {self._filetype}"
                )
        else:
            logging.info(f"{self._filetype} not found")
        logging.info(f"Download {self._filetype}")
        data = self.fetch_url(self._data_url)
        logging.debug(f"Write {self._filetype} to {self._data}")
        self._save_data(data)
        return self._load_data()


class Metadata_updater(Data_updater):
    """Update local metadata for local dataset storage for BUSCO markers.

    This writes to user defined cfg_dir for storing metata about version and
    dataset names for downloading BUSCO markers.
    """

    def __init__(self, database_url):
        """Initialize metadata object and local db file."""
        super().__init__(database_url)
        self._filetype = "metadata"
        self._data = phyling.config.cfg_dir / "metadata.pickle"
        self._data_url = f"{database_url}/file_versions.tsv"
        self._metadata_md5_url = f"{database_url}/file_versions.tsv.hash"

    @property
    def get_metadata_path(self):
        """Obtain the local metadata DB path."""
        return self._data

    @property
    def _get_local_md5(self):
        """Retrieve the local cached md5 checksums."""
        with open(self._data, "rb") as f:
            _, md5 = pickle.load(f)
        return md5

    @property
    def _get_remote_md5(self):
        """Retrieve the remote md5 checksum for a file."""
        return self.fetch_url(self._metadata_md5_url).decode().strip("\n")

    def _load_data(self):
        """Load local metadata DB file."""
        with open(self._data, "rb") as f:
            data, _ = pickle.load(f)
        return data

    def _save_data(self, data):
        """Save metadata into local DB file."""
        markerset = {}
        for line in data.decode().split("\n"):
            line = line.split("\t")
            if line[-1] == "lineages":
                markerset[line[0]] = {
                    "url": f"{self._database_url}/lineages/{line[0]}.{line[1]}.tar.gz",
                    "md5": line[2],
                }
        with open(self._data, "wb") as f:
            pickle.dump((markerset, hashlib.md5(data).hexdigest()), f)


class HMM_markerset_updater(Data_updater):
    """Update and store local copy of HMM files into destination."""

    def __init__(self, database_url, output_dir, metadata: dict, name: str):
        """Initialize the markerset objects for download and retrieval."""
        super().__init__(database_url)
        self._filetype = "HMM markerset"
        self._database_url = "/".join([self._database_url, "lineages"])
        self._markerset = metadata[name]
        self._data = output_dir / name / "md5sum"  # Using md5sum file to represent the database.
        self._data_url = self._markerset["url"]

    @property
    def _get_local_md5(self):
        """Read and retrieve local md5 checksums."""
        with open(self._data) as f:
            md5 = f.read().strip("\n")
        return md5

    @property
    def _get_remote_md5(self):  # Actually is the md5sum recorded in markerset_dict
        """Get the md5sum stored in theh markerset_dict object."""
        return self._markerset["md5"]

    def _save_data(self, data):
        """Save the marker dataset into a folder for HMM use."""
        output = self._data.parent
        output.mkdir(parents=True, exist_ok=True)
        logging.debug(f"Save to {output}.tar.gz")
        with open(f"{output}.tar.gz", "wb") as f:
            f.write(data)
        logging.debug(f"Gunzip the {self._filetype} to {output}")
        with tarfile.open(f"{output}.tar.gz", "r:gz") as f:
            f.extractall(output.parent)

        md5 = hashlib.md5(open(f"{output}.tar.gz", "rb").read()).hexdigest()
        logging.debug(f"Write the md5 checksum to {output}/md5sum")
        print(md5, file=open(self._data, "w"))
        logging.debug(f"Remove {output}.tar.gz")
        Path(f"{output}.tar.gz").unlink()


# def md5(file):
#     with open(file, 'rb') as f:
#         hasher = hashlib.md5()
#         chunk = f.read(8195)
#     while chunk:
#         hasher.update(chunk)
#         chunk = f.read(8195)
#     return hasher.hexdigest()


def download(markerset, **kwargs) -> None:
    """
    Help to download/update BUSCO v5 markerset to a local folder.

    First it checks whether the metadata file is exist under the config folder ~/.phyling. A missing or outdated file
    will trigger the module to download/update the metadata.

    Passing "list" to markerset argument will list all the available/already downloaded markersets. Passing a valid
    name to the markerset argument will download the markerset to the config folder ~/.phyling/HMM.
    """
    metadata_updater = Metadata_updater(database_url=phyling.config.database)
    markerset_dict = metadata_updater.updater()

    if markerset == "list":
        # Get the dictionary with database as key and convert it into a list
        url_list = [hmm_markerset for hmm_markerset in markerset_dict.keys()]
        url_list.sort()
        exist_markerset = []
        for markerset in Path(phyling.config.cfg_dir, phyling.config.default_HMM).iterdir():
            if markerset.name in url_list and markerset.is_dir():
                exist_markerset.append(markerset.name)
        exist_markerset.sort()
        print("Available datasets:\n")

        # Adjust databases display according to the terminal size
        width, _ = shutil.get_terminal_size((80, 24))
        col = width // 40
        url_list = [url_list[x : x + col] for x in range(0, len(url_list), col)]
        col_width = max(len(word) for row in url_list for word in row) + 3  # padding
        for row in url_list:
            # Print the database list
            print(" ".join(word.ljust(col_width) for word in row))

        print()
        print("Downloaded datasets:\n")
        exist_markerset = [exist_markerset[x : x + col] for x in range(0, len(exist_markerset), col)]
        for row in exist_markerset:
            # Print the database list
            print(" ".join(word.ljust(col_width) for word in row))
    else:
        hmm_markerset_updater = HMM_markerset_updater(
            database_url=phyling.config.database,
            output_dir=Path(phyling.config.cfg_dir, phyling.config.default_HMM),
            metadata=markerset_dict,
            name=markerset,
        )
        hmm_markerset_updater.updater()
