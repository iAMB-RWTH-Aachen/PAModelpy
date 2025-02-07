from Scripts.enzymatic_data_analysis import create_membrane_enz_complex, find_var_for_complex
import pandas as pd
import pytest
def test_create_membrane_enz_complex():
    # Arrange
    memprot_dict = {
                        'P0AGM7': [{'CE_URAt2pp_P0AGM7': {'f': 22.0, 'b': 0}}, 14],
                        'P75892': [{'CE_URAt2pp_P75892': {'f': 22.0, 'b': 0}}, 12]
                    }
    expected_list = ['P0AGM7', 'P75892']

    # Act
    actual_list = create_membrane_enz_complex(memprot_dict)

    # Assert
    assert actual_list == expected_list

@pytest.mark.parametrize(
    "enz_complex, enzyme_db, p_conc_expected, alpha_expected",
    [
        (
            'P0AGM6_P75895',
            pd.DataFrame({
                'uniprotID': ['P0AGM7', 'P75892'],
                'SP4+TX100_avg': [20, 87],
                'alpha_helix_units': [27, 30]
            }),
            0,
            0
        ),
        # Add more test cases here
        (
            'P0AGM7_P75892',
            pd.DataFrame({
                'uniprotID': ['P0AGM7', 'P75892'],
                'SP4+TX100_avg': [30, 70],
                'alpha_helix_units': [28, 32]
            }),
            30,
            32
        ),
        # Add as many test cases as needed
    ]
)
def test_find_var_for_complex(enz_complex, enzyme_db, p_conc_expected, alpha_expected):
    # Act
    p_conc_actual, alpha_actual = find_var_for_complex(enz_complex, enzyme_db)

    # Assert
    assert p_conc_actual == p_conc_expected
    assert alpha_actual == alpha_expected