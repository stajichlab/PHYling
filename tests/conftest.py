import shutil
import tarfile

import pytest

import phyling.config as config


@pytest.fixture(scope="module")
def pack_current_metadata():
    metadata = config.cfg_dir / config.metadata
    data_folder = config.cfg_dir / config.default_HMM
    temp_tar = config.cfg_dir / "temp.tgz"
    with tarfile.open(temp_tar, "w") as tar:
        if metadata.exists() and metadata.is_file():
            tar.add(metadata, arcname=metadata.name)
            metadata.unlink()
        if data_folder.exists() and data_folder.is_dir():
            tar.add(data_folder, arcname=data_folder.name)
            shutil.rmtree(data_folder)
    yield
    if data_folder.exists() and data_folder.is_dir():
        shutil.rmtree(data_folder)
    metadata.unlink(missing_ok=True)
    with tarfile.open(temp_tar, "r") as tar:
        tar.extractall(config.cfg_dir, filter="fully_trusted")
    temp_tar.unlink()


@pytest.fixture(scope="class")
def shared_tmpdir_class(tmp_path_factory: pytest.TempPathFactory):
    return tmp_path_factory.mktemp("shared_tmpdir")
