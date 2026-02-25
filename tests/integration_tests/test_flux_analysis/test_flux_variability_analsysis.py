import pandas as pd
import pytest
from cobra import Reaction

from src.PAModelpy.flux_analysis import flux_variability_analysis
from src.PAModelpy import Enzyme, EnzymeVariable
from tests.unit_tests.test_pamodel.test_pamodel import build_toy_pam



@pytest.fixture
def toy_pam():
    pam = build_toy_pam(sensitivity = False)
    pam.change_reaction_bounds('R1', lower_bound=-1e3, upper_bound=1e3)
    return pam


@pytest.mark.parametrize("variable_type, attr_name", [
    (Reaction, 'reactions'),
    (Enzyme, 'enzyme_variables'),
    (EnzymeVariable, 'enzyme_variables')
]
                         )
def test_fva_runs_for_different_object_types(toy_pam, variable_type, attr_name):
    df = flux_variability_analysis(model = toy_pam,
                                   variable_type = variable_type
                                   )
    assert len(df) == len(getattr(toy_pam, attr_name))

def test_fva_gets_same_result_in_parallel_as_without(toy_pam):
    df_ref = flux_variability_analysis(model = toy_pam
                                   )
    df = flux_variability_analysis(model = toy_pam,
                                   processes = 2
                                   )
    for i, row in df.iterrows():
        assert i in df_ref.index
        assert df_ref['minimum'].loc[i] == pytest.approx(row['minimum'], abs=1e-6)
        assert df_ref['maximum'].loc[i] == pytest.approx(row['maximum'], abs=1e-6)