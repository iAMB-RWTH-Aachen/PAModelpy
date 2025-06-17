import pytest
from cobra import Reaction

from tests.unit_tests.tests_pamodel import build_toy_pam
from src.PAModelpy.utils.recombinant_protein_overexpression import *
from src.PAModelpy import Enzyme

@pytest.fixture
def mock_pam():
    return build_toy_pam()

def test_add_protein_export_to_pam(mock_pam):
    rxn = add_protein_export(mock_pam, protein_name="testprot")
    assert rxn.id == "testprot_production"
    assert "testprot_c" in [m.id for m in rxn.metabolites]
    assert any(r.id == "testprott" for r in mock_pam.reactions)
    assert any(b.id == "EX_testprot_e" for b in mock_pam.reactions)

def test_add_aminoacid_sequence_to_production_rxns(mock_pam):
    reaction = Reaction("fake_rxn")
    mock_pam.add_reactions([reaction])
    aa_freq = {'met__L_c': 0.5, 'val__L_c': 0.5}
    add_aminoacid_sequence_to_production_rxns(mock_pam, aa_freq, reaction)
    assert len(reaction.metabolites) == 2
    for m in reaction.metabolites:
        assert m.id in mock_pam.metabolites

def test_add_recombinant_protein_production_and_export(dummy_sequence_file, mock_pam):
    rxn = add_recombinant_protein_production_and_export(aa_txt_file=dummy_sequence_file,
                                                        pam=mock_pam,
                                                        protein_name="recombinase",
                                                        export_efficiency=0.9,
                                                        molecular_weight=50000)
    assert rxn.id == "recombinase_production"
    assert "recombinase" in mock_pam.enzyme_variables
    assert "recombinase_c" in [m.id for m in rxn.metabolites]
    assert any(r.id == "recombinase_production" for r in mock_pam.reactions)
