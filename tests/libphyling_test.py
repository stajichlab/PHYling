from pathlib import Path

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyling.libphyling import bp_mrtrans, msa_generator, trim_gaps


class TestMSAGenerator:
    test_inputs = [x for x in (Path("example") / "pep").iterdir()]
    inputs_with_duplicates = [Path(x) for x in ["Species_A.fasta.gz", "Species_A.fasta", "Species_B.fasta.gz"]]

    def test_msa_generator(self):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            msa_generator(self.inputs_with_duplicates)
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1

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

    def test_filter_orthologs_with_insufficient_orthologs(self, msa_instance):
        # Simulate having small orthologs in the instance
        msa_instance.orthologs = {
            "hmm_1": {"species_1_gene_100", "species_2_gene_211"},
            "hmm_2": {"species_1_gene_105", "species_2_gene_223"},
            "hmm_3": {"species_1_gene_155", "species_2_gene_253"},
            "hmm_4": {"species_3_gene_344", "species_4_gene_466"},
        }
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            msa_instance.filter_orthologs()
        assert len(msa_instance.orthologs) == 0  # No orthologs should remain
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 1

    def test_filter_orthologs_attribute_error(self, msa_instance):
        # Simulate not having orthologs attribute in the instance
        with pytest.raises(AttributeError):
            msa_instance.filter_orthologs()


def create_msa(records, ids):
    return MultipleSeqAlignment([SeqRecord(Seq(rec), id=id) for rec, id in zip(records, ids)])


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
    pep_msa = create_msa(
        [
            "MI-T-C",
            "M-AT-C",
            "MI-TG-",
            "M-GT-C",
        ],
        [
            "Species_A",
            "Species_B",
            "Species_C",
            "Species_D",
        ],
    )
    cds_msa = create_msa(
        [
            "ATGATC---ACG---TGT",
            "ATG---GCTACT---TGT",
            "ATGATA---ACAGGA---",
            "ATG---GGTACT---TGC",
        ],
        [
            "Species_A",
            "Species_B",
            "Species_C",
            "Species_D",
        ],
    )

    def test_trim_gaps_with_peptide(self):
        result_msa = trim_gaps(self.pep_msa, gaps=0.6)
        assert len(result_msa) == 4
        assert str(result_msa[0].seq) == "MI-TC"
        assert str(result_msa[1].seq) == "M-ATC"
        assert str(result_msa[2].seq) == "MI-T-"
        assert str(result_msa[3].seq) == "M-GTC"

    def test_trim_gaps_with_peptide_no_trim(self):
        result_msa = trim_gaps(self.pep_msa, gaps=0.8)
        assert len(result_msa) == 4
        assert str(result_msa[0].seq) == "MI-T-C"
        assert str(result_msa[1].seq) == "M-AT-C"
        assert str(result_msa[2].seq) == "MI-TG-"
        assert str(result_msa[3].seq) == "M-GT-C"

    def test_trim_gaps_with_cds_msa(self):
        result_msa = trim_gaps(self.pep_msa, self.cds_msa, gaps=0.6)
        assert len(result_msa) == 4
        assert str(result_msa[0].seq) == "ATGATC---ACGTGT"
        assert str(result_msa[1].seq) == "ATG---GCTACTTGT"
        assert str(result_msa[2].seq) == "ATGATA---ACA---"
        assert str(result_msa[3].seq) == "ATG---GGTACTTGC"

    def test_trim_gaps_with_cds_msa_no_trim(self):
        result_msa = trim_gaps(self.pep_msa, self.cds_msa, gaps=0.8)
        assert len(result_msa) == 4
        assert str(result_msa[0].seq) == "ATGATC---ACG---TGT"
        assert str(result_msa[1].seq) == "ATG---GCTACT---TGT"
        assert str(result_msa[2].seq) == "ATGATA---ACAGGA---"
        assert str(result_msa[3].seq) == "ATG---GGTACT---TGC"

    def test_trim_gaps_at_threshold_boundary(self):
        result_msa = trim_gaps(self.pep_msa, gaps=0.5)
        assert len(result_msa) == 4
        assert str(result_msa[0].seq) == "MTC"
        assert str(result_msa[1].seq) == "MTC"
        assert str(result_msa[2].seq) == "MT-"
        assert str(result_msa[3].seq) == "MTC"

    def test_trim_gaps_invalid_gaps_value(self):
        with pytest.raises(ValueError):
            trim_gaps(self.pep_msa, gaps=1.5)

    def test_trim_gaps_msa_metadata(self):
        result_msa = trim_gaps(self.pep_msa, gaps=0.5)
        assert len(result_msa) == 4
        assert str(result_msa[0].id) == "Species_A"
        assert str(result_msa[1].id) == "Species_B"
        assert str(result_msa[2].id) == "Species_C"
        assert str(result_msa[3].id) == "Species_D"
