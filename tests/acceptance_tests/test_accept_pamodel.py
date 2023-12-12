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

from Scripts.pam_generation import set_up_toy_pam, set_up_ecoli_pam

def test_toy_model_sensitivity_coefficients_relations_sum_is_correct():
    #arrange
    toy_pam = set_up_toy_pam()
    #act
    toy_pam.test(-0.001)
    #assert
    assert_sensitivity_coefficients(toy_pam)

def test_ecoli_model_sensitivity_coefficients_relations_sum_is_correct():
    # arrange
    ecoli_pam = set_up_ecoli_pam()
    # act
    ecoli_pam.test()
    # assert
    assert_sensitivity_coefficients(ecoli_pam)


#######################################################################################################
#HELPER METHODS
#######################################################################################################

def assert_sensitivity_coefficients(pamodel):
    feasibility_deviation = pamodel.solver.problem.Params.FeasibilityTol * len(pamodel.variables) *4

    # calculate fraction of protein space occupied
    sum_enzymes = pamodel.calculate_sum_of_enzymes()
    phi_e0 = pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].ub
    alpha = sum_enzymes / phi_e0

    # calculate sum of coefficients
    sum_FCSC = calculate_sum_of_flux_capacity_sensitivity_coefficients(pamodel)
    sum_ECSC = calculate_sum_of_enzyme_capacity_sensitivity_coefficients(pamodel)
    PCSC = calculate_proteome_capacity_sensitivity_coefficients(pamodel)
    sum_ESC = calculate_sum_of_enzyme_sensitivity_coefficients(pamodel)

    # validation sums
    enzyme_flux_sum = sum_ESC + sum_FCSC + (1 - alpha) * PCSC
    enzyme_sum = -sum_ESC + sum_ECSC + alpha * PCSC

    #assert
    assert 1 == pytest.approx(sum_ECSC+sum_FCSC + PCSC, feasibility_deviation)
    # assert 1 == pytest.approx(enzyme_flux_sum, feasibility_deviation)
    # assert 0 == pytest.approx(enzyme_sum, feasibility_deviation)

def calculate_sum_of_enzymes(pamodel):
    sum = 0  # mg/gcdw/h
    for enzyme in pamodel.enzyme_variables:
        sum += enzyme.concentration
    return sum

def calculate_sum_of_flux_capacity_sensitivity_coefficients(pamodel):
    return sum(pamodel.capacity_sensitivity_coefficients[
                       (pamodel.capacity_sensitivity_coefficients['constraint'] == 'flux_ub') | (
                                   pamodel.capacity_sensitivity_coefficients['constraint'] == 'flux_lb')].coefficient)

def calculate_sum_of_enzyme_capacity_sensitivity_coefficients(pamodel):
    return sum(pamodel.capacity_sensitivity_coefficients[
                       (pamodel.capacity_sensitivity_coefficients['constraint'] == 'enzyme_min') | (
                                   pamodel.capacity_sensitivity_coefficients[
                                       'constraint'] == 'enzyme_max')].coefficient)

def calculate_proteome_capacity_sensitivity_coefficients(pamodel):
    return pamodel.capacity_sensitivity_coefficients[
        pamodel.capacity_sensitivity_coefficients['constraint'] == 'proteome'].coefficient.iloc[0]

def calculate_sum_of_enzyme_sensitivity_coefficients(pamodel):
    return sum(pamodel.enzyme_sensitivity_coefficients.coefficient)

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

def build_toy_pam():
    Config.BIOMASS_REACTION = 'R7'
    Config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    Config.CO2_EXHANGE_RXNID = 'R8'
    Config.ACETATE_EXCRETION_RXNID = 'R9'

    model = load_json_model('toy_model.json')
    active_enzyme = build_active_enzyme_sector(Config)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme,
                      unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config)
    return pamodel