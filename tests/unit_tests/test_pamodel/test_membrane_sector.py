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
from src.PAModelpy.MembraneSector import MembraneSector, UnusedMembraneSector

def test_unused_membrane_sector_updates_variable_side_correctly():
    # Arrange
    sut = build_toy_model(membrane_sector=True, unused_membrane_sector=False)

    membrane_constraint = sut.constraints["membrane"]

    ups_sector = sut.sectors.get_by_id("UnusedEnzymeSector")
    lin_rxn = sut.reactions.get_by_id(ups_sector.id_list[0])

    forward_var = lin_rxn.forward_variable
    reverse_var = lin_rxn.reverse_variable

    # Baseline (before)
    baseline_coeffs = membrane_constraint.get_linear_coefficients(
        [forward_var, reverse_var]
    )

    baseline_fwd = baseline_coeffs[forward_var]
    baseline_rev = baseline_coeffs[reverse_var]

    # Act
    unused_membrane_sector = UnusedMembraneSector(sut)

    # Re-read coefficients (after)
    updated_coeffs = membrane_constraint.get_linear_coefficients(
        [forward_var, reverse_var]
    )

    updated_fwd = updated_coeffs[forward_var]
    updated_rev = updated_coeffs[reverse_var]

    # Expected delta behavior
    expected_slope = unused_membrane_sector.slope

    expected_fwd = baseline_fwd + expected_slope
    expected_rev = baseline_rev - expected_slope

    # Assert full transformation, not just slope correctness
    assert updated_fwd == pytest.approx(expected_fwd)
    assert updated_rev == pytest.approx(expected_rev)

def test_unused_membrane_sector_updates_membrane_constraint_upper_bound_correctly():
    # Arrange
    sut = build_toy_model(membrane_sector=True, unused_membrane_sector=False) # Build a model with membrane sector without unused membrane sector
    original_ub = sut.constraints['membrane'].ub

    # Act
    unused_membrane_sector = UnusedMembraneSector(sut) # Add unused membrane sector which automatically links the unused membrane sector to the original membrane sector
    new_ub = sut.constraints['membrane'].ub
    
    assert new_ub < original_ub # The new intercept will always be smaller than the original intercept because a portion of the available S/V ratio is allocated for unused membrane enyzmes (substracted). 
    assert new_ub == pytest.approx(
        original_ub - unused_membrane_sector.intercept
    )

@pytest.mark.parametrize('membrane_complex', ['E1_E2_E3', 'E4', 'E5_E6', 'E7_E8_E9_E10'])
def test_if_membrane_enzyme_is_excluded_when_updating_membrane_sector(membrane_complex):
    # Arrange
    enzyme_complex = membrane_complex

    # Act
    sut = build_toy_model(membrane_sector=True)
    enzyme_variable = sut.enzyme_variables[0]
    sut.constraints[sut.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients( # Added a dummy membrane enzyme to the tpc constraint
        {  
                enzyme_variable.forward_variable: 1,
                enzyme_variable.reverse_variable: 1,
        }
    )
    sut.sectors.get_by_id('MembraneSector')._update_membrane_constraint(0.02, sut)

    # Assert
    tpc = sut.constraints[sut.TOTAL_PROTEIN_CONSTRAINT_ID]
    var_names = {v.name for v in tpc.expression.free_symbols}
    membrane_sector = sut.sectors.get_by_id("MembraneSector")

    if enzyme_complex in membrane_sector.membrane_proteins:
        assert not any(
            name.startswith(enzyme_complex)
            for name in var_names
        )

def test_if_get_alpha_number_for_enz_complex_works():
    # Arrange
    membrane_sector = build_membrane_sector()
    enz_complex = "E7_E8_E9_E10"
    expected_alpha_number = 192

    # Act
    actual_alpha_number = membrane_sector._get_alpha_number_for_enz_complex(enz_complex)

    # Assert
    assert expected_alpha_number == actual_alpha_number

def test_if_add_membrane_constraint_works():
    # Arrange
    toy_pam = build_toy_model(membrane_sector=True)
    alpha_numbers_dict = {'E1_E2_E3': 20,
                          'E7_E8_E9_E10': 48,
                          'E14': 12}

    ##calling out the variable names
    biomass_var = toy_pam.reactions.get_by_id(toy_pam.BIOMASS_REACTION).forward_variable

    for enz_complex in toy_pam.enzyme_variables:
        if enz_complex.id in alpha_numbers_dict.keys():
            expected_coefficients = {
                biomass_var: (-toy_pam.membrane_sector.slope),
                enz_complex.forward_variable: (1e-6 * alpha_numbers_dict[enz_complex.id]
                               * toy_pam.membrane_sector.area_alpha
                               * toy_pam.membrane_sector.unit_factor),
                enz_complex.reverse_variable: (1e-6 * alpha_numbers_dict[enz_complex.id]
                               * toy_pam.membrane_sector.area_alpha
                               * toy_pam.membrane_sector.unit_factor)
            }

    # Act
    ## Retrieve the membrane constraint
    membrane_constraint = toy_pam.constraints['membrane']
    print(f"Membrane constraint: {membrane_constraint}")
    for enz_complex in toy_pam.enzyme_variables:
        retrieved_coefficients = membrane_constraint.get_linear_coefficients([biomass_var, enz_complex.forward_variable, enz_complex.reverse_variable])

    #Debug
    print('Expected coefficients: ', expected_coefficients)
    print('Retrieved coefficients: ', retrieved_coefficients)

    # Assert
    assert membrane_constraint is not None, "Membrane constraint not found"
    assert retrieved_coefficients == expected_coefficients, f"Expected {expected_coefficients}, but got {retrieved_coefficients}"
    assert membrane_constraint.lb == 0
    assert membrane_constraint.ub == toy_pam.membrane_sector.intercept

##helper methods
def build_membrane_sector(unused_membrane_sector: bool=False):
    alpha_numbers_dict = {"E1": 20,
                          "E2": 20,
                          "E3": 20,
                          "E4": 8,
                          "E5": 1,
                          "E6": 1,
                          "E7": 48,
                          "E8": 48,
                          "E9": 48,
                          "E10": 48,
                          "E11": 12,
                          "E12": 12,
                          "E13": 12,
                          "E14": 12}

    enzyme_location = {"E1": "Cell inner membrane",
                       "E2": "Cell inner membrane",
                       "E3": "Cell inner membrane",
                       "E4": "Cytosol",
                       "E5": "Cytosol",
                       "E6": "Cytosol",
                       "E7": "Cell inner membrane",
                       "E8": "Cell inner membrane",
                       "E9": "Cell inner membrane",
                       "E10": "Cell inner membrane",
                       "E11": "Cytosol",
                       "E12": "Cytosol",
                       "E13": "Cytosol",
                       "E14": "Cell inner membrane"}
    # mu = 0.1, 0=0.005
    membrane_sector = MembraneSector(sv_slope=-0.1, sv_0=1,
                                         alpha_numbers_dict=alpha_numbers_dict,
                                         enzyme_location=enzyme_location, usable_area_fraction=1,
                                         unused_membrane_sector=unused_membrane_sector)

    return membrane_sector

def build_toy_model(sensitivity:bool=True, membrane_sector: bool=False, unused_membrane_sector: bool =False):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'R11'
    config.GLUCOSE_EXCHANGE_RXNID = 'R2'
    config.CO2_EXHANGE_RXNID = 'R12'
    config.ACETATE_EXCRETION_RXNID = 'R7'

    nmbr_reactions = 12

    # Building Active Enzyme Sector
    kcat_fwd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 10e-2, 10e-2]
    kcat_rev = [kcat for kcat in kcat_fwd]
    enz_complex_id = ['E1_E2_E3', 'E4', 'E5_E6', 'E7_E8_E9_E10', 'E11', 'E12', 'E13', 'E14']
    gpr_string = [['g1','g2','g3'], ['g4'],['g5','g6'],['g7','g8','g9','g10'],['g11','g12'],['g13'],['g14'],['g15']]
    pra_string = ['g1 and g2 and g3', 'g4', 'g5 and g6', 'g7 and g8 and g9 and g10', 'g11 or g12', 'g13', 'g14', 'g15']

    rxn2protein = {}
    protein2gene = {}
    for i in range(nmbr_reactions-4): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+2}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass
        rxn2protein = {**rxn2protein, **{rxn_id: {f'{enz_complex_id[i]}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6),
                                                             'molmass': 1e6, 'genes': gpr_string[i],
                                                             'protein_reaction_association': pra_string[i]}}}}
        protein2gene = {**protein2gene, **{f'E{i+2}': gpr_string[i]}}

    active_enzyme = ActiveEnzymeSector(rxn2protein = rxn2protein, protein2gene=protein2gene, configuration=config)

    # Building Translational Protein Sector
    translation_enzyme = TransEnzymeSector(id_list = ['R11'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=config)

    # Building Unused Enzyme Sector
    unused_enzyme = UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=config)

    # Building Membrane Sector
    if membrane_sector:
        membrane_sector = build_membrane_sector(unused_membrane_sector=unused_membrane_sector)
    else:
        membrane_sector = None

    # Building the toy_pam
    model = load_json_model('Models/toy_model.json')
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints',
                      sensitivity=sensitivity,
                      active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme,
                      membrane_sector=membrane_sector,
                      p_tot=0.2, configuration=config)

    return pamodel

def assert_bounds(model_ori, model_copy):
    for key, var in model_ori.variables.items():
        assert var.ub == model_copy.variables[key].ub
        assert var.lb == model_copy.variables[key].lb

    for key, cons in model_ori.constraints.items():
        assert cons.ub == model_copy.constraints[key].ub
        assert cons.lb == model_copy.constraints[key].lb

def assert_total_protein_content(model_ori, model_copy):
    assert model_ori.p_tot == model_copy.p_tot
    tot_prot_cons_id = model_ori.TOTAL_PROTEIN_CONSTRAINT_ID
    assert model_ori.constraints[tot_prot_cons_id].ub == model_copy.constraints[tot_prot_cons_id].ub