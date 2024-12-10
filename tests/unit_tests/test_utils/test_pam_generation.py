import os.path

import cobra.io
import numpy as np
import pandas as pd
import pytest

from Scripts.toy_ec_pam import build_toy_gem
from src.utils.pam_generation import parse_reaction2protein, set_up_pam

def test_if_rxn2protein_info_is_correctly_parsed():
    # Arrange
    toy_enzyme_db = pd.DataFrame({'rxn_id':['R1', 'R2', 'R2', 'R3', 'R2'],
                                  'enzyme_id': ['E1', 'E2a', 'E2b_E2c', 'E3', 'E2a'],
                                  'gene': [['gene1'],['gene2a'], ['gene2b', 'gene2c'], ['gene3'], ['gene2a']],
                                  'GPR': ['gene1','gene2a or (gene2b and gene2c)', 'gene2a or (gene2b and gene2c)',
                                          'gene3', 'gene2a or (gene2b and gene2c)'],
                                  'molMass':[1, np.nan,2,3,2],
                                  'kcat_values':[1,2,np.nan,3,1.5],
                                  'direction':['f','f', 'f', 'f', 'b']
                                  }
                                 )
    toy_model = build_toy_gem()

    expected_rxn2protein = {
        'R1':  {
            'E1': {
                'enzyme_id': 'E1', 'f': 1.0, 'b': 0,
                'genes': ['gene1'], 'protein_reaction_association': [['E1']],
                'molmass': 1.0}
        },
        'R2': {
            'E2a': {
                'enzyme_id': 'E2a', 'f': 2.0, 'b': 1.5,
                'genes': ['gene2a'], 'protein_reaction_association': [['E2a']],
                'molmass': 39959.4825},
            'E2b_E2c': {
                'enzyme_id': 'E2b_E2c', 'f': 0, 'b': 0,
                'genes': ['gene2b', 'gene2c'], 'protein_reaction_association': [['E2b', 'E2c']],
                'molmass': 2.0}
        },
        'R3': {
            'E3': {
                'enzyme_id': 'E3', 'f': 3.0, 'b': 0,
                'genes': ['gene3'], 'protein_reaction_association': [['E3']],
                'molmass': 3.0}
        }
    }
    expected_protein2gpr = {'E1': [['gene1']], 'E2a': [['gene2a']], 'E2b_E2c': [['gene2b', 'gene2c']], 'E3': [['gene3']]}

    # Apply
    rxn2protein, protein2gpr = parse_reaction2protein(toy_enzyme_db, toy_model)

    # Assert
    for output_dict, expected_dict in zip([rxn2protein, protein2gpr], [expected_rxn2protein, expected_protein2gpr]):
        assert all([expected_dict[key] == value for key, value in output_dict.items()])

def test_if_set_up_pam_can_build_ecolicore_pam():
    #Arrange
    pam_data_file = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_core.xlsx')
    ecolicore_gem = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))

    #Apply
    ecolicore_pam = set_up_pam(pam_data_file,
                               ecolicore_gem,
                               total_protein = 0.1699,
                               sensitivity=False,
                               adjust_reaction_ids=True)

    ecolicore_pam.optimize()

    #Assert
    assert ecolicore_pam.objective.value > 0

def test_if_set_up_pam_can_build_iML1515():
    #Arrange
    pam_data_file = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_241209.xlsx')
    iml1515 = os.path.join('Models', 'iML1515.xml')

    #Apply
    pam = set_up_pam(pam_data_file,
                               iml1515,
                               sensitivity=False,
                               adjust_reaction_ids=True)

    pam.optimize()

    #Assert
    assert pam.objective.value > 0