from phyling.libphyling import dict_merge


def test_dict_merge():
    input = [
        {
            "hmm_1": "species_1_gene_100", "hmm_2": "species_1_gene_105", "hmm_3": "species_1_gene_155",
        },
        {
            "hmm_1": "species_2_gene_211", "hmm_2": "species_2_gene_223", "hmm_3": "species_2_gene_253",
        },
        {
            "hmm_1": "species_3_gene_311", "hmm_2": "species_3_gene_323", "hmm_4": "species_3_gene_344",
        },
        {
            "hmm_1": "species_4_gene_414", "hmm_3": "species_4_gene_432", "hmm_4": "species_4_gene_466",
        },
    ]
    expect = {
        "hmm_1": set(["species_1_gene_100", "species_2_gene_211", "species_3_gene_311", "species_4_gene_414"]),
        "hmm_2": set(["species_1_gene_105", "species_2_gene_223", "species_3_gene_323"]),
        "hmm_3": set(["species_1_gene_155", "species_2_gene_253", "species_4_gene_432"]),
        "hmm_4": set(["species_3_gene_344", "species_4_gene_466"]),
    }
    assert dict_merge(input) == expect
