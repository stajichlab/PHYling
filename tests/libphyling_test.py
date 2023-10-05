from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyling.libphyling import bp_mrtrans, dict_merge, msa_generator, trim_gaps


class TestMSAGenerator:
    test_inputs = [x for x in (Path("example") / "pep").iterdir()]

    @pytest.fixture
    def msa_instance(self):
        return msa_generator(self.test_inputs)

    def test_filter_orthologs_with_orthologs(self, msa_instance):
        # Simulate having orthologs in the instance
        msa_instance.orthologs = {
            "hmm_1": {"species_1_gene_100", "species_2_gene_211", "species_3_gene_311", "species_4_gene_414"},
            "hmm_2": {"species_1_gene_105", "species_2_gene_223", "species_3_gene_323"},
            "hmm_3": {"species_1_gene_155", "species_2_gene_253", "species_4_gene_432"},
            "hmm_4": {"species_3_gene_344", "species_4_gene_466"},
        }
        msa_instance.filter_orthologs()
        assert len(msa_instance.orthologs) == 3
        assert "hmm_4" not in msa_instance.orthologs
        assert len(msa_instance.orthologs["hmm_1"]) == 4  # Number of samples

    def test_filter_orthologs_with_small_orthologs(self, msa_instance):
        # Simulate having small orthologs in the instance
        msa_instance.orthologs = {
            "hmm_1": {"species_1_gene_100", "species_2_gene_211"},
            "hmm_2": {"species_1_gene_105", "species_2_gene_223"},
            "hmm_3": {"species_1_gene_155", "species_2_gene_253"},
            "hmm_4": {"species_3_gene_344", "species_4_gene_466"},
        }
        msa_instance.filter_orthologs()
        assert len(msa_instance.orthologs) == 0  # No orthologs should remain

    def test_filter_orthologs_attribute_error(self, msa_instance):
        # Simulate not having orthologs attribute in the instance
        with pytest.raises(AttributeError):
            msa_instance.filter_orthologs()


def create_msa(records):
    return MultipleSeqAlignment([SeqRecord(Seq(rec)) for rec in records])


class Testdictmerge:
    def test_dict_merge_empty_input(self):
        result = dict_merge([])
        assert result == {}

    def test_dict_merge_single_dict(self):
        dicts_list = [{"a": 1, "b": 2}]
        result = dict_merge(dicts_list)
        assert result == {"a": {1}, "b": {2}}

    def test_dict_merge_multiple_dicts(self):
        dicts_list = [{"a": 1, "b": 2}, {"a": 3, "c": 4}]
        result = dict_merge(dicts_list)
        assert result == {"a": {1, 3}, "b": {2}, "c": {4}}

    def test_dict_merge_duplicate_values(self):
        dicts_list = [{"a": 1, "b": 2}, {"a": 1, "c": 3}]
        result = dict_merge(dicts_list)
        assert result == {"a": {1}, "b": {2}, "c": {3}}

    def test_dict_merge_empty_dicts(self):
        dicts_list = [{}, {}]
        result = dict_merge(dicts_list)
        assert result == {}


class Testbpmrtrans:
    pep_msa = create_msa(["-MSLR-L-", "-M-LRQL-", "-MS--QL-"])

    def test_bp_mrtrans_basic(self):
        cds_seqs = [
            SeqRecord(Seq("ATGTCATTGCGACTA")),
            SeqRecord(Seq("ATGTTGCGACAACTA")),
            SeqRecord(Seq("ATGTCACAACTA")),
        ]
        results = bp_mrtrans(pep_msa=self.pep_msa, cds_seqs=cds_seqs)
        assert str(results[0].seq) == "---ATGTCATTGCGA---CTA---"
        assert str(results[1].seq) == "---ATG---TTGCGACAACTA---"
        assert str(results[2].seq) == "---ATGTCA------CAACTA---"

    def test_bp_mrtrans_with_stop_codon(self):
        cds_seqs = [
            SeqRecord(Seq("ATGTCATTGCGACTA")),
            SeqRecord(Seq("ATGTGATTGCGACAACTA")),
            SeqRecord(Seq("ATGTCACAACTA")),
        ]
        results = bp_mrtrans(pep_msa=self.pep_msa, cds_seqs=cds_seqs)
        assert str(results[0].seq) == "---ATGTCATTGCGA---CTA---"
        assert str(results[1].seq) == "---ATG---TTGCGACAACTA---"
        assert str(results[2].seq) == "---ATGTCA------CAACTA---"


from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyling.libphyling import bp_mrtrans, dict_merge, msa_generator, trim_gaps


class TestMSAGenerator:
    test_inputs = [x for x in (Path("example") / "pep").iterdir()]

    @pytest.fixture
    def msa_instance(self):
        return msa_generator(self.test_inputs)

    def test_filter_orthologs_with_orthologs(self, msa_instance):
        # Simulate having orthologs in the instance
        msa_instance.orthologs = {
            "hmm_1": {"species_1_gene_100", "species_2_gene_211", "species_3_gene_311", "species_4_gene_414"},
            "hmm_2": {"species_1_gene_105", "species_2_gene_223", "species_3_gene_323"},
            "hmm_3": {"species_1_gene_155", "species_2_gene_253", "species_4_gene_432"},
            "hmm_4": {"species_3_gene_344", "species_4_gene_466"},
        }
        msa_instance.filter_orthologs()
        assert len(msa_instance.orthologs) == 3
        assert "hmm_4" not in msa_instance.orthologs
        assert len(msa_instance.orthologs["hmm_1"]) == 4  # Number of samples

    def test_filter_orthologs_with_small_orthologs(self, msa_instance):
        # Simulate having small orthologs in the instance
        msa_instance.orthologs = {
            "hmm_1": {"species_1_gene_100", "species_2_gene_211"},
            "hmm_2": {"species_1_gene_105", "species_2_gene_223"},
            "hmm_3": {"species_1_gene_155", "species_2_gene_253"},
            "hmm_4": {"species_3_gene_344", "species_4_gene_466"},
        }
        msa_instance.filter_orthologs()
        assert len(msa_instance.orthologs) == 0  # No orthologs should remain

    def test_filter_orthologs_attribute_error(self, msa_instance):
        # Simulate not having orthologs attribute in the instance
        with pytest.raises(AttributeError):
            msa_instance.filter_orthologs()


def create_msa(records, ids):
    return MultipleSeqAlignment([SeqRecord(Seq(rec), id=id) for rec, id in zip(records, ids)])


class Testdictmerge:
    def test_dict_merge_empty_input(self):
        result = dict_merge([])
        assert result == {}

    def test_dict_merge_single_dict(self):
        dicts_list = [{"a": 1, "b": 2}]
        result = dict_merge(dicts_list)
        assert result == {"a": {1}, "b": {2}}

    def test_dict_merge_multiple_dicts(self):
        dicts_list = [{"a": 1, "b": 2}, {"a": 3, "c": 4}]
        result = dict_merge(dicts_list)
        assert result == {"a": {1, 3}, "b": {2}, "c": {4}}

    def test_dict_merge_duplicate_values(self):
        dicts_list = [{"a": 1, "b": 2}, {"a": 1, "c": 3}]
        result = dict_merge(dicts_list)
        assert result == {"a": {1}, "b": {2}, "c": {3}}

    def test_dict_merge_empty_dicts(self):
        dicts_list = [{}, {}]
        result = dict_merge(dicts_list)
        assert result == {}


class Testbpmrtrans:
    pep_msa = create_msa(["-MSLR-L-", "-M-LRQL-", "-MS--QL-"], ["Species_A", "Species_B", "Species_C"])

    def test_bp_mrtrans_basic(self):
        cds_seqs = [
            SeqRecord(Seq("ATGTCATTGCGACTA"), id="Species_A"),
            SeqRecord(Seq("ATGTTGCGACAACTA"), id="Species_B"),
            SeqRecord(Seq("ATGTCACAACTA"), id="Species_C"),
        ]
        results = bp_mrtrans(pep_msa=self.pep_msa, cds_seqs=cds_seqs)
        assert str(results[0].seq) == "---ATGTCATTGCGA---CTA---"
        assert str(results[1].seq) == "---ATG---TTGCGACAACTA---"
        assert str(results[2].seq) == "---ATGTCA------CAACTA---"

    def test_bp_mrtrans_with_stop_codon(self):
        cds_seqs = [
            SeqRecord(Seq("ATGTCATTGCGACTA")),
            SeqRecord(Seq("ATGTGATTGCGACAACTA")),
            SeqRecord(Seq("ATGTCACAACTA")),
        ]
        results = bp_mrtrans(pep_msa=self.pep_msa, cds_seqs=cds_seqs)
        assert str(results[0].seq) == "---ATGTCATTGCGA---CTA---"
        assert str(results[1].seq) == "---ATG---TTGCGACAACTA---"
        assert str(results[2].seq) == "---ATGTCA------CAACTA---"

    def test_bp_mrtrans_with_cds_msa_id(self):
        cds_seqs = [
            SeqRecord(Seq("ATGTCATTGCGACTA"), id="Species_A_cds"),
            SeqRecord(Seq("ATGTTGCGACAACTA"), id="Species_B_cds"),
            SeqRecord(Seq("ATGTCACAACTA"), id="Species_C_cds"),
        ]
        results = bp_mrtrans(pep_msa=self.pep_msa, cds_seqs=cds_seqs)
        assert str(results[0].id) == "Species_A_cds"
        assert str(results[1].id) == "Species_B_cds"
        assert str(results[2].id) == "Species_C_cds"

    # def test_bp_mrtrans_wo_cds_msa_id(self):
    #     cds_seqs = [
    #         SeqRecord(Seq("ATGTCATTGCGACTA")),
    #         SeqRecord(Seq("ATGTTGCGACAACTA")),
    #         SeqRecord(Seq("ATGTCACAACTA")),
    #     ]
    #     results = bp_mrtrans(pep_msa=self.pep_msa, cds_seqs=cds_seqs)
    #     assert str(results[0].id) == "Species_A"
    #     assert str(results[1].id) == "Species_B"
    #     assert str(results[2].id) == "Species_C"


class Testtrimgaps:
    pep_msa1 = create_msa(["-MG--A", "M-GT-C", "MMGTG-"], ["Species_A", "Species_B", "Species_C"])
    pep_msa2 = create_msa(["MTG--A", "M-GT-C"], ["Species_A", "Species_B"])

    def test_trim_gaps_basic(self):
        result_msa = trim_gaps(self.pep_msa1, gaps=0.5)
        assert len(result_msa) == 3
        assert str(result_msa[0].seq) == "-MG-A"
        assert str(result_msa[1].seq) == "M-GTC"
        assert str(result_msa[2].seq) == "MMGT-"

    def test_trim_gaps_with_cds_msa(self):
        cds_msa = create_msa(["ATGACTGGA------GCT", "ATG---GCTACT---TGT"], ["Species_A", "Species_B"])
        result_msa = trim_gaps(self.pep_msa2, cds_msa, gaps=0.5)
        assert len(result_msa) == 2
        assert str(result_msa[0].seq) == "ATGGGAGCT"
        assert str(result_msa[1].seq) == "ATGGCTTGT"

    def test_trim_gaps_at_threshold_boundary(self):
        result_msa = trim_gaps(self.pep_msa2, gaps=0.5)
        assert len(result_msa) == 2
        assert str(result_msa[0].seq) == "MGA"
        assert str(result_msa[1].seq) == "MGC"

    def test_trim_gaps_no_trim(self):
        result_msa = trim_gaps(self.pep_msa2, gaps=0.51)
        assert len(result_msa) == 2
        assert str(result_msa[0].seq) == "MTG-A"
        assert str(result_msa[1].seq) == "M-GTC"

    def test_trim_gaps_invalid_gaps_value(self):
        with pytest.raises(ValueError):
            trim_gaps(self.pep_msa2, gaps=1.5)

    def test_trim_gaps_msa_metadata(self):
        result_msa = trim_gaps(self.pep_msa1, gaps=0.5)
        assert len(result_msa) == 3
        assert str(result_msa[0].id) == "Species_A"
        assert str(result_msa[1].id) == "Species_B"
        assert str(result_msa[2].id) == "Species_C"
