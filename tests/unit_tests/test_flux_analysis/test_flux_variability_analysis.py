import pandas as pd
import pytest
from cobra.core import Reaction
from tests.unit_tests.test_tamodel.test_tamodel import build_toy_tam

from src.flux_analysis import flux_variability_analysis
from src.PAModelpy import Enzyme
from src.TAModelpy import Transcript

def test_flux_variablity_analysis_for_reactions():
    # Arrange
    sut = build_toy_tam()
    sut.change_reaction_bounds(rxn_id='R1',
                               lower_bound=0, upper_bound=1e-1)
    reactions_to_test = ['R2', 'R3', 'R4']
    expected_results = pd.DataFrame({'minimum': [0.079742,0,0.026581], 'maximum': [0.079742,0,0.026581]},
                                    index=reactions_to_test)

    # Act
    fva_results = flux_variability_analysis(model = sut,
                                            variable_type = Reaction,
                                            variable_list = reactions_to_test)

    # Assert
    assert pd.testing.assert_frame_equal(expected_results, fva_results, check_exact=False, atol=1e-3) is None

def test_flux_variablity_analysis_for_enzymes():
    # Arrange
    sut = build_toy_tam()
    sut.change_reaction_bounds(rxn_id='R1',
                               lower_bound=0, upper_bound=1e-1)
    enzymes_to_test = ['E2', 'E3', 'E4']
    expected_results = pd.DataFrame({'minimum': [0.159484,0,0.026581], 'maximum': [0.159484,0,0.026581]},
                                    index=enzymes_to_test)

    # Act
    fva_results = flux_variability_analysis(model = sut,
                                            variable_type = Enzyme,
                                            variable_list = enzymes_to_test)

    # Assert
    assert pd.testing.assert_frame_equal(expected_results, fva_results, check_exact=False, atol=1e-3) is None

def test_flux_variablity_analysis_for_transcripts():
    # Arrange
    sut = build_toy_tam()
    sut.change_reaction_bounds(rxn_id='R1',
                               lower_bound=0, upper_bound=1e-1)
    transcripts_to_test = ['mRNA2', 'mRNA3', 'mRNA4']
    expected_results = pd.DataFrame({'minimum': [7.655235e-11,0,1.275872e-11], 'maximum': [3.987101e-09,0,6.645169e-10]},
                                    index=transcripts_to_test)
    # Act
    fva_results = flux_variability_analysis(model = sut,
                                            variable_type = Transcript,
                                            variable_list = transcripts_to_test)
    # Assert
    assert pd.testing.assert_frame_equal(expected_results, fva_results, check_exact=False, atol=1e-3) is None
