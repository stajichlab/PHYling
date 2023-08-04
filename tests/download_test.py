from phyling.download import Metadata_updater, HMM_markerset_updater

database = "https://busco-data.ezlab.org/v5/data"


def test_download_module(tmpdir):
    metadata_updater = Metadata_updater(database_url=database, cfg_dir=tmpdir)
    markerset_dict = metadata_updater.updater()

    hmm_markerset_updater = HMM_markerset_updater(
        database_url=database,
        output_dir="tests/HMM",
        metadata=markerset_dict,
        name="fungi_odb10",
    )
    hmm_markerset_updater.updater()
