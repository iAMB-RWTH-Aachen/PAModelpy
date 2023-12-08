import unittest
import pytest
from cobra.io import load_json_model
import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(
            os.path.dirname( #testing dir
                os.path.dirname(__file__))))) #this dir
print(os.getcwd())

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
    print(sut.constraints[constraint_name + 'b'].get_linear_coefficients([rxn.reverse_variable]))
    sut.change_kcat_value(enzyme_id, kcats ={'R1':{'f': input_kcat, 'b': input_kcat}})
    print('works until here')
    print(sut.constraints[constraint_name+'b'].get_linear_coefficients([rxn.reverse_variable]))
    coeff_b = sut.constraints[constraint_name+'b'].get_linear_coefficients([rxn.reverse_variable])[rxn.reverse_variable]

    #/(3600*1e-6) to correct for dimension modifications in the model
    model_kcat_b = 1/coeff_b/(3600*1e-6)
    coeff_f = sut.constraints[constraint_name+'f'].get_linear_coefficients([rxn.forward_variable])[rxn.forward_variable]
    model_kcat_f = 1/coeff_f/(3600*1e-6)

    #assert
    assert input_kcat == pytest.approx(model_kcat_b, 1e-4)
    assert input_kcat == pytest.approx(model_kcat_f, 1e-4)


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
