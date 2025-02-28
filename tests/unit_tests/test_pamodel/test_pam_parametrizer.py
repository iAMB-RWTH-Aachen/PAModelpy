import unittest
import pytest
from cobra.io import load_json_model
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.dirname(
            os.path.dirname( #testing dir
                os.path.dirname(__file__))))) #this dir

from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.configuration import Config
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from Scripts.create_pamodel_from_diagnostics_file import create_pamodel_from_diagnostics_file
from src.PAModelpy.utils.pam_generation import _order_enzyme_complex_id

@pytest.mark.parametrize('case, enzyme_id, expected_kcat_f',
                         [
                             ('case_01',
                              'E1_E2_E3',
                              1004),

                             ('case_02',
                              'E1_E2_E3',
                              1004),

                             ('case_02',
                              'E10_E7_E8_E9',
                              1014),

                              ('case_02',
                              'E10_E7_E8_E9',
                              1014)
                         ])
def test_if_create_pamodel_from_diagnostics_file_is_correct(case, enzyme_id, expected_kcat_f):
    # Arrange
    toy_pam = build_toy_model()

    # Act
    data_path = 'tests/data/diagnostics_file_for_toy.xlsx'
    toy_pam = create_pamodel_from_diagnostics_file(data_path, toy_pam, case)
    kcats = toy_pam.enzyme_variables.get_by_id(enzyme_id).kcats.values()
    actual_kcat_f = list(kcats)[0]['f'] # Getting only the forward kcat, because only f kcat is in the diagnostics file

    # Assert
    assert actual_kcat_f == expected_kcat_f

### Helper functions ###

def build_toy_model(sensitivity:bool=True, membrane_sector: bool=False):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'R11'
    config.GLUCOSE_EXCHANGE_RXNID = 'R2'
    config.CO2_EXHANGE_RXNID = 'R12'
    config.ACETATE_EXCRETION_RXNID = 'R7'

    nmbr_reactions = 12

    # Building Active Enzyme Sector
    kcat_fwd = [1000, 1000, 2000, 1500, 2500, 1122, 1023, 1055, 650]
    kcat_rev = [kcat for kcat in kcat_fwd]
    enz_complex_id = ['E1_E2_E3', 'E4', 'E5_E6', 'E7_E8_E9_E10', 'E11', 'E12', 'E13', 'E14']
    gpr_string = [['g1','g2','g3'], ['g4'],['g5','g6'],['g7','g8','g9','g10'],['g11','g12'],['g13'],['g14'],['g15']]
    pra_string = ['g1 and g2 and g3', 'g4', 'g5 and g6', 'g7 and g8 and g9 and g10', 'g11 or g12', 'g13', 'g14', 'g15']

    rxn2protein = {}
    protein2gene = {}

    i=0
    for enz_complex in enz_complex_id:
        enz_id = _order_enzyme_complex_id(enz_complex)
        enz_complex_id[i] = enz_id
        i+=1
    print(enz_complex_id)

    for i in range(nmbr_reactions-4): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+2}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass
        rxn2protein = {**rxn2protein, **{rxn_id: {f'{enz_complex_id[i]}':{'enyzme_id': enz_complex_id[i], 'f': kcat_fwd[i], 'b': kcat_rev[i],
                                                             'molmass': 1e6, 'genes': gpr_string[i],
                                                             'protein_reaction_association': pra_string[i]}}}}
        protein2gene = {**protein2gene, **{f'E{i+2}': gpr_string[i]}}
    active_enzyme = ActiveEnzymeSector(rxn2protein = rxn2protein, protein2gene=protein2gene, configuration=config)

    # Building Translational Protein Sector
    translation_enzyme = TransEnzymeSector(id_list = ['R11'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=config)

    # Building Unused Enzyme Sector
    unused_enzyme = UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=config)

    # Building the toy_pam
    model = load_json_model('Models/toy_model.json')
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints',
                      sensitivity=sensitivity,
                      active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme,
                      p_tot=0.2, configuration=config)

    return pamodel


if __name__ == '__main__':
    toy_pam = build_toy_model()
    toy_pam.objective = 'R11'
    toy_pam.optimize()
    kcats = toy_pam.enzyme_variables.get_by_id('E1_E2_E3').kcats.values()
    value_f = list(kcats)[0]['f']
    print(value_f)