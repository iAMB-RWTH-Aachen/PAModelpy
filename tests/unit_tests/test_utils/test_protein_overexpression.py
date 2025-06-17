import pytest
import tempfile

from src.PAModelpy.utils.recombinant_protein_overexpression import *

@pytest.fixture
def dummy_sequence_file():
    aa_seq = "MKTFFVVAL"  # A short 1-letter amino acid sequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
        f.write(aa_seq)
        f_path = f.name
    yield f_path

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
