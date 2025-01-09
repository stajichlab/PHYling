from __future__ import annotations

import shutil
import tarfile

import pytest

import phyling.libphyling
from phyling.pipeline.download import download


class PackMetadata:
    def __init__(self) -> None:
        self.metadata = phyling.libphyling.METADATA_FILE
        self.data_folder = phyling.libphyling.HMM_DIR
        self.temp_tar = phyling.CFG_DIR / "temp.tgz"

    def pack_metadata(self):
        with tarfile.open(self.temp_tar, "w") as tar:
            if self.metadata.exists() and self.metadata.is_file():
                tar.add(self.metadata, arcname=self.metadata.name)
                self.metadata.unlink()
            if self.data_folder.exists() and self.data_folder.is_dir():
                tar.add(self.data_folder, arcname=self.data_folder.name)
                shutil.rmtree(self.data_folder)

    def unpack_metadata(self):
        if self.data_folder.exists() and self.data_folder.is_dir():
            shutil.rmtree(self.data_folder)
        self.metadata.unlink(missing_ok=True)
        with tarfile.open(self.temp_tar, "r") as tar:
            tar.extractall(phyling.CFG_DIR, filter="fully_trusted")
        self.temp_tar.unlink()


@pytest.fixture(scope="class")
def hide_metadata_wo_download():
    obj = PackMetadata()
    obj.pack_metadata()
    yield
    obj.unpack_metadata()


@pytest.fixture(scope="module")
def hide_metadata_with_download():
    obj = PackMetadata()
    obj.pack_metadata()
    download("poxviridae_odb10")
    yield
    obj.unpack_metadata()


@pytest.fixture(scope="class")
def shared_tmpdir_class(tmp_path_factory: pytest.TempPathFactory):
    return tmp_path_factory.mktemp("shared_tmpdir")


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="run slow tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
