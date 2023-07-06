import logging
import pickle
import shutil
import sys
import hashlib
import tarfile
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError


class metadata_updater():
    """
    A metadata_updater. Will check the metadata pickle file and download the metadata if it is absent/outdated.

    Attributes
    ----------
    database : str
        The root url of the BUSCO database.
    metadata_pickle : str
        The path of the pickle file. (Default: "file_versions.pickle")

    Methods
    -------
    check_metadata()
        Check the metadata pickle file to see if download/update is needed.
    """
    def __init__(self, cfg_dir: Path, database: str, metadata_pickle="file_versions.pickle"):
        self.__database = database
        self.__metadata_md5_url = f'{database}/file_versions.tsv.hash'
        self.__metadata_pickle = cfg_dir / metadata_pickle

    def check_metadata(self) -> dict:
        """
        Check the metadata pickle file to see if download/update is needed.
        
        Return
        ------
        dict: {"markerset"}
        """
        logging.debug("Check matadata pickle file exist or not ...")
        if self.__metadata_pickle.is_file():
            logging.debug("Metadata pickle file exist.")
            with open(self.__metadata_pickle, 'rb') as f:
                metadata = pickle.load(f)
            
            logging.debug("Compare md5 checksum between online version and the md5 recorded in the pickle file...")
            if fetch_url(self.__metadata_md5_url).decode().strip('\n') == metadata['md5']:
                logging.debug("md5 is the same. No need to update metadata")
                return metadata['markerset']
            else:
                logging.info("The local file md5 is different from the online version")
        else:
            logging.info("Metadata not found")
        logging.info("Update metadata")
        metadata = self.__build_metadata()
        with open(self.__metadata_pickle, 'wb') as f:
            pickle.dump(metadata, f)
        return metadata['markerset']

    def __build_metadata(self) -> dict:
        """
        Download the metadata (file_versions.tsv) from BUSCO and build a dictionary which contains the md5 checksum of the 
        metadata itself and a dictionary with available BUSCO markersets. The dictionary is also saved as a pickle file
        which provides a quick check for the next execution.

        Return
        ------
        dict: {"md5", "markerset"}
        """
        metadata = fetch_url(f'{self.__database}/file_versions.tsv')
        md5 = hashlib.md5(metadata).hexdigest()
        logging.debug(f'Generating markerset dictionary ...')
        markerset = {}
        for line in metadata.decode().split('\n'):
            line = line.split('\t')
            if line[-1] == 'lineages':
                markerset[line[0]] = {"url": f'{self.__database}/lineages/{".".join([line[0], line[1], "tar.gz"])}',
                                      "md5": line[2]}
        return {"md5": md5, "markerset": markerset}

def fetch_url(url) -> bytes:
    """
    Download the content from the url.

    Attributes
    ----------
    url : str
        The url source.

    Return
    ------
    bytes
    """
    try:
        logging.debug(f'Download {url} ...')
        with urlopen(url) as response:
            content = response.read()
    except HTTPError:
        logging.error("URL not found or currently unavailble")
        sys.exit(1)
    return content

# def md5(file):
#     with open(file, 'rb') as f:
#         hasher = hashlib.md5()
#         chunk = f.read(8195)
#     while chunk:
#         hasher.update(chunk)
#         chunk = f.read(8195)
#     return hasher.hexdigest()

def download(args) -> None:
    # Check whether the file_version.pickle is exist. A missing/outdated file will trigger downloading of the file_version.tsv
    # and convert to a pickle file with a md5 checksum. Will return a dictionary with database and its md5 checksum at the end.
    obj_metadata_updater = metadata_updater(cfg_dir=args.cfg_dir, database=args.database)
    markerset_dict = obj_metadata_updater.check_metadata()

    if args.markerset == 'list':
        # Get the dictionary with database as key and convert it into a list
        url_list = [hmm_markerset for hmm_markerset in markerset_dict.keys()]
        url_list.sort()
        print('Available datasets:\n')

        # Adjust databases display according to the terminal size
        width, _ = shutil.get_terminal_size((80, 24))
        col = width // 40
        url_list = [url_list[x: x + col] for x in range(0, len(url_list), col)]
        col_width = max(len(word) for row in url_list for word in row) + 3  # padding
        for row in url_list:
            # Print the database list
            print(" ".join(word.ljust(col_width) for word in row))
        
    else:
        logging.debug("Check hmm markerset exist or not ...")
        output = args.output / args.markerset
        output.mkdir(parents=True, exist_ok=True)
        if Path(f'{output}/md5sum').is_file():
            logging.debug("Hmm markerset exist")
            with open(f'{output}/md5sum', 'r') as f:
                md5 = f.read().strip('\n')
            if md5 == markerset_dict[args.markerset]['md5']:
                logging.info("md5 is the same. No need to update hmm markerset")
                sys.exit(0)
            else:
                logging.info("The local file md5 is different from the online version")
                logging.info("Remove the local file")
                shutil.rmtree(f'{output}')

        logging.info("Update hmm markerset")
        markerset_content = fetch_url(markerset_dict[args.markerset]['url'])
        logging.debug(f'Write to {output}.tar.gz')
        with open(f'{output}.tar.gz', 'wb') as f:
            f.write(markerset_content)
        logging.debug(f'Gunzip the hmm markerset to {output}')
        with tarfile.open(f'{output}.tar.gz', 'r:gz') as f:
            f.extractall()

        md5 = hashlib.md5(open(f'{output}.tar.gz', 'rb').read()).hexdigest()
        logging.debug(f'Write the md5 checksum to {output}/md5sum')
        print(md5, file=open(f'{output}/md5sum', 'w'))
        Path(f'{output}.tar.gz').unlink()
        logging.info(f'Hmm markerset is saved in {output}')
