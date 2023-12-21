import unittest
import pytest
from cobra.io import load_json_model
import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(
            os.path.dirname( #testing dir
                os.path.dirname(__file__))))) #this dir

from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.configuration import Config
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector



def test_if_pamodel_change_kcat_function_works():
    #arrange
    sut = build_toy_pam()
    input_kcat = 10
    enzyme_id = 'E1'
    rxn= sut.reactions.get_by_id('R1')
    constraint_name = 'EC_E1_'

    #act
    sut.change_kcat_value(enzyme_id, kcats ={'R1':{'f': input_kcat, 'b': input_kcat}})
    coeff_b = sut.constraints[constraint_name+'b'].get_linear_coefficients([rxn.reverse_variable])[rxn.reverse_variable]
    #/(3600*1e-6) to correct for dimension modifications in the model
    model_kcat_b = 1/coeff_b/(3600*1e-6)

    coeff_f = sut.constraints[constraint_name+'f'].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
    # /(3600*1e-6) to correct for dimension modifications in the model
    model_kcat_f = 1/coeff_f/(3600*1e-6)

    #assert
    assert input_kcat == pytest.approx(model_kcat_b, 1e-4)
    assert input_kcat == pytest.approx(model_kcat_f, 1e-4)

def test_if_pamodel_change_reaction_bounds_function_works_with_sensitivity_constraints():
    #arrange
    toy_pam = build_toy_pam(sensitivity=True)
    new_R1_ub = 500
    new_R1_lb = 50

    #act
    toy_pam.change_reaction_bounds(rxn_id='R1', lower_bound= new_R1_lb, upper_bound=new_R1_ub)
    new_toy_R1_lb = -toy_pam.constraints['R1_lb'].ub
    new_toy_R1_ub = toy_pam.constraints['R1_ub'].ub

    #assert
    assert new_R1_lb == new_toy_R1_lb
    assert new_R1_ub == new_toy_R1_ub


def test_if_pamodel_change_reaction_bounds_function_works_without_sensitivity_constraints():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    new_R1_ub = 500
    new_R1_lb = 50

    # act
    toy_pam.change_reaction_bounds(rxn_id='R1', lower_bound=new_R1_lb, upper_bound=new_R1_ub)
    new_toy_R1_lb = toy_pam.reactions.get_by_id('R1').lower_bound
    new_toy_R1_ub = toy_pam.reactions.get_by_id('R1').upper_bound

    # assert
    assert new_R1_lb == new_toy_R1_lb
    assert new_R1_ub == new_toy_R1_ub


def test_if_pamodel_change_total_protein_function_works():
    #arrange
    toy_pam = build_toy_pam()
    protein_fraction_sectors = toy_pam.p_tot*1e3 - toy_pam.constraints[toy_pam.TOTAL_PROTEIN_CONSTRAINT_ID].ub
    new_ptot = 0.5

    #act
    toy_pam.change_total_protein_constraint(new_ptot)
    new_toy_ptot = toy_pam.constraints[toy_pam.TOTAL_PROTEIN_CONSTRAINT_ID].ub

    #assert
    assert new_ptot*1e3-protein_fraction_sectors == new_toy_ptot

def test_if_pamodel_change_enzyme_bounds_function_with_sensitivity_works():
    # arrange
    toy_pam = build_toy_pam(sensitivity=True)
    new_E1_max = 0.5
    new_E1_min = 0.005

    # act
    toy_pam.change_enzyme_bounds(enzyme_id='E1', lower_bound=new_E1_min, upper_bound=new_E1_max)
    new_toy_E1_min = -toy_pam.constraints['E1_min'].ub
    new_toy_E1_max = toy_pam.constraints['E1_max'].ub

    # assert
    assert new_E1_min == new_toy_E1_min
    assert new_E1_max == new_toy_E1_max

def test_if_pamodel_change_enzyme_bounds_function_without_sensitivity_works():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    new_E1_max = 0.5
    new_E1_min = 0.005
    # act
    toy_pam.change_enzyme_bounds(enzyme_id='E1', lower_bound=new_E1_min, upper_bound=new_E1_max)
    new_toy_E1_min = toy_pam.enzyme_variables.get_by_id('E1').lower_bound
    new_toy_E1_max = toy_pam.enzyme_variables.get_by_id('E1').upper_bound

    # assert
    assert new_E1_min == new_toy_E1_min
    assert new_E1_max == new_toy_E1_max

def test_if_pamodel_copy_function_works():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    # act
    toy_pam_copy = toy_pam.copy()
    # assert
    assert_bounds(toy_pam, toy_pam_copy)
    assert_total_protein_content(toy_pam, toy_pam_copy)

def test_if_pamodel_enzyme_get_kcat_values_correctly():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    enzyme_to_return = 'E3'
    reaction_to_return = 'R3'
    kcat_to_return = 1/(3600*1e-6) #unit correction to be consistent with the way it is saved (see toymodel setup function)

    # act
    kcat_returned = toy_pam.enzymes.get_by_id(enzyme_to_return).get_kcat_values([reaction_to_return])

    # assert
    assert kcat_to_return == kcat_returned['f']


#######################################################################################################
#HELPER METHODS
#######################################################################################################
def build_active_enzyme_sector(Config):
    n=9
    kcat_fwd = [1, 0.5, 1, 1, 0.5 ,0.45, 1.5]  # High turnover for exchange reactions
    kcat_rev = [kcat for kcat in kcat_fwd]
    rxn2kcat = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id: {f'E{i+1}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6), 'molmass': 1e6}}}}

    return ActiveEnzymeSector(rxn2protein = rxn2kcat, configuration=Config)

def build_unused_protein_sector(Config):
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=Config)

def build_translational_protein_sector(Config):
    return TransEnzymeSector(id_list = ['R7'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=Config)

def build_toy_pam(sensitivity = True):
    Config.BIOMASS_REACTION = 'R7'
    Config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    Config.CO2_EXHANGE_RXNID = 'R8'
    Config.ACETATE_EXCRETION_RXNID = 'R9'

    model = load_json_model('tests/unit_tests/toy_model.json')
    active_enzyme = build_active_enzyme_sector(Config)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme, sensitivity = sensitivity,
                      unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config)
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