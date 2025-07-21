import pytest
import tempfile


from Scripts.pam_generation_uniprot_id import set_up_ecoli_pam
from src.PAModelpy.utils.recombinant_protein_expression import (_get_substrate_uptake_rate_for_fixed_growth_rate,
                                                                read_sequence_from_file,
                                                                match_aminoacid_to_model_identifiers_and_frequency
                                                                )

@pytest.fixture
def dummy_sequence_file():
    aa_seq = "MKTFFVVAL"  # A short 1-letter amino acid sequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
        f.write(aa_seq)
        f_path = f.name
    yield f_path

@pytest.fixture
def mock_pam():
    return set_up_ecoli_pam(sensitivity=False)

def test_read_sequence_from_file(dummy_sequence_file):
    seq = read_sequence_from_file(dummy_sequence_file)
    assert isinstance(seq, str)
    assert seq == "MKTFFVVAL"


def test_match_aminoacid_to_model_identifiers_and_frequency():
    seq = "MKTFFVVAL"
    aa_freq = match_aminoacid_to_model_identifiers_and_frequency(seq)

    expected_keys = {"met__L_c", "lys__L_c", "thr__L_c", "phe__L_c", "val__L_c", "ala__L_c", "leu__L_c"}
    assert set(aa_freq.keys()) == expected_keys
    assert abs(sum(aa_freq.values()) - 1.0) < 1e-6  # should sum to 1


def test_get_subtrate_uptake_rate_for_fixed_growth_rate(mock_pam):
    substrate_rxn = mock_pam.reactions.EX_glc__D_e
    uptake_id = substrate_rxn.id

    rate = _get_substrate_uptake_rate_for_fixed_growth_rate(
        pam=mock_pam,
        substrate_uptake_id=uptake_id,
        growth_rate=0.1
    )

    assert rate < 0  # uptake should be negative
    assert isinstance(rate, float)

