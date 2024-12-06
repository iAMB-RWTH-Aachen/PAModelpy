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
from src.PAModelpy.MembraneSector import MembraneSector

def test_if_add_membrane_constraint_works():
    # Arrange
    toy_pam = build_toy_pam(membrane_sector=build_membrane_sector())

        #calling out the variable names
    biomass_var = toy_pam.reactions.get_by_id(toy_pam.BIOMASS_REACTION).forward_variable
    forward_var = toy_pam.enzyme_variables.P43660_P43661_P43662.forward_variable
    reverse_var = toy_pam.enzyme_variables.P43660_P43661_P43662.reverse_variable

    expected_coefficients = {
        biomass_var: (-toy_pam.membrane_sector.slope),
        forward_var: (1e-6 * toy_pam.membrane_sector.alpha_numbers_dict['P43661']
                       * toy_pam.membrane_sector.area_alpha
                       * toy_pam.membrane_sector.unit_factor),
        reverse_var: (1e-6 * toy_pam.membrane_sector.alpha_numbers_dict['P43661']
                       * toy_pam.membrane_sector.area_alpha
                       * toy_pam.membrane_sector.unit_factor)
    }

    # Act
    toy_pam.add_membrane_constraint()

    # Retrieve the membrane constraint
    membrane_constraint = toy_pam.constraints['membrane']
    print(f"Membrane constraint: {membrane_constraint}")
    retrieved_coefficients = membrane_constraint.get_linear_coefficients([biomass_var, forward_var, reverse_var])

    #Debug
    print('Expected coefficients: ', expected_coefficients)
    print('Retrieved coefficients: ', retrieved_coefficients)

    # Assert
    assert membrane_constraint is not None, "Membrane constraint not found"
    assert retrieved_coefficients == expected_coefficients, f"Expected {expected_coefficients}, but got {retrieved_coefficients}"
    assert membrane_constraint.lb == 0
    assert membrane_constraint.ub == toy_pam.membrane_sector.intercept

def test_if_membrane_constraints_the_same_for_sensitivity_with_and_without_toypam():

    # building the mcpam models
    membrane_sector = build_membrane_sector()
    sut_toypam_true = build_toy_pam(sensitivity=True, membrane_sector=membrane_sector)
    sut_toypam_false = build_toy_pam(sensitivity=False, membrane_sector=membrane_sector)

    # calling out the membrane constraints for both mcpam models
    membrane_constraint_sut_toypam_true = sut_toypam_true.constraints['membrane']
    membrane_constraint_sut_toypam_false = sut_toypam_false.constraints['membrane']

    # calling out the upper and lower bound of the constraints for both mcpam models
    models = [sut_toypam_true, sut_toypam_false]
    df_true = pd.DataFrame(columns=['lb', 'ub'])
    df_false = pd.DataFrame(columns=['lb', 'ub'])

    for model in models:
        lb = []
        ub = []
        if model.sensitivity:
            for constraint in model.constraints:
                lb.append(constraint.lb)
                ub.append(constraint.ub)

            df_true['lb'] = lb
            df_true['ub'] = ub

        else:
            for constraint in model.constraints:
                lb.append(constraint.lb)
                ub.append(constraint.ub)

            df_false['lb'] = lb
            df_false['ub'] = ub

    assert membrane_constraint_sut_toypam_true == membrane_constraint_sut_toypam_false
    assert df_true['lb'] == df_false['lb']
    assert df_true['ub'] == df_false['ub']

##helper methods
def build_active_enzyme_sector(Config):
    n=9
    kcat_fwd = [1, 0.5, 1, 1, 0.5 ,0.45, 1.5]  # High turnover for exchange reactions
    kcat_rev = [kcat for kcat in kcat_fwd]
    gpr_string = [['gene1', 'gene2', 'gene3'], ['gene2', 'gene3'], ['gene3'], ['gene4'], ['gene5'], ['gene6'], ['gene7'], ['gene8']]
    pra_string = [[['P43660', 'P43661', 'P43662']], [['P43661', 'P43662']], [['P43662']], [['P43663']], [['P43664']], [['P43665']], [['P43666']]]
    rxn2protein = {}
    protein2gene = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2protein = {**rxn2protein,
                       **{rxn_id: {f'P4366{i}': {'f': kcat_fwd[i] / (3600 * 1e-6), 'b': kcat_rev[i] / (3600 * 1e-6),
                                                 'molmass': 1e6, 'genes': gpr_string[i],
                                                 'protein_reaction_association': pra_string[i]}}}}
        protein2gene = {**protein2gene, **{f'P4366{i}': gpr_string[i]}}

    return ActiveEnzymeSector(rxn2protein = rxn2protein, protein2gene=protein2gene, configuration=Config)

def build_unused_protein_sector(Config):
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=Config)

def build_translational_protein_sector(Config):
    return TransEnzymeSector(id_list = ['R7'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=Config)

def build_membrane_sector():
    alpha_numbers_dict = {"P43660": 12,
                          "P43661": 8,
                          "P43662": 1,
                          "P43663": 12,
                          "P43664": 40,
                          "P43665": 15,
                          "P43666": 15}

    enzyme_location = {"P43660": "Unknown",
                       "P43661": "Cell membrane",
                       "P43662": "Cell membrane",
                       "P43663": "Cytoplasm",
                       "P43664": "Cytoplasm",
                       "P43665": "Cell membrane",
                       "P43666": "Cell membrane"}

    return MembraneSector(area_avail_mu=6.2129, area_avail_0=4.7522,
                                     alpha_numbers_dict=alpha_numbers_dict, enzyme_location=enzyme_location)

def build_toy_pam(sensitivity = True, membrane_sector = None):
    Config.BIOMASS_REACTION = 'R7'
    Config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    Config.CO2_EXHANGE_RXNID = 'R8'
    Config.ACETATE_EXCRETION_RXNID = 'R9'

    model = load_json_model('Data/toy_model_cell_membrane_localization.json')
    active_enzyme = build_active_enzyme_sector(Config)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme, sensitivity = sensitivity,
                      unused_sector = unused_enzyme, membrane_sector=membrane_sector,
                      p_tot=0.6*1e-3, configuration=Config)
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