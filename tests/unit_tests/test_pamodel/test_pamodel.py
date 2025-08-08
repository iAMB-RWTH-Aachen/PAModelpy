import pytest
from cobra.io import load_json_model

from src.PAModelpy import PAModel,Config,ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector, CatalyticEvent
from Scripts.pam_generation_uniprot_id import set_up_ecoli_pam, set_up_ecolicore_pam
from tests.unit_tests.test_pamodel.test_pam_generation_ecoli import set_up_toy_pam_with_isozymes_and_enzymecomplex


def test_if_pamodel_change_kcat_function_works():
    #arrange
    sut = build_toy_pam(sensitivity=False)
    input_kcat = 10
    enzyme_id = 'E1'
    rxn= sut.reactions.get_by_id('CE_R1_'+enzyme_id)
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

def test_if_pamodel_gets_reaction_id_from_complex_catalytic_reaction_id():
    # Arrange
    sut = build_toy_pam(sensitivity=False)
    ce_reaction_id = "CE_CYTBO3_4pp_P0ABJ1_P0ABJ5_P0ABJ7_P0ABJ8"
    reaction_id = "CYTBO3_4pp"

    # Act
    parsed_reaction_id = CatalyticEvent._extract_reaction_id_from_catalytic_reaction_id(ce_reaction_id)

    # Assert
    assert parsed_reaction_id == reaction_id

def test_if_pamodel_change_kcat_function_works_with_catalytic_reactions():
    #arrange
    sut = set_up_ecoli_pam(sensitivity=False)
    input_kcat = 10
    enzyme_id = 'P0ABJ1'
    rxn_id = "CYTBO3_4pp"
    ce_rxn= sut.reactions.query(f'CE_{rxn_id}_{enzyme_id}')[0]
    enzyme_complex_id = "_".join(ce_rxn.id.split("_")[3:])
    constraint_name = f'EC_{enzyme_complex_id}_'

    #act
    sut.change_kcat_value(enzyme_complex_id, kcats ={ce_rxn.id:{'f': input_kcat, 'b': input_kcat}})
    coeff_b = sut.constraints[constraint_name+'b'].get_linear_coefficients([ce_rxn.reverse_variable])[ce_rxn.reverse_variable]
    #/(3600*1e-6) to correct for dimension modifications in the model
    model_kcat_b = 1/coeff_b/(3600*1e-6)

    coeff_f = sut.constraints[constraint_name+'f'].get_linear_coefficients([ce_rxn.forward_variable])[ce_rxn.forward_variable]
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

def test_if_pamodel_change_enzyme_bounds_function_can_change_when_lb_geq_ub():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    E1_min, E1_max = 0.005, 0.5
    toy_pam.change_enzyme_bounds(enzyme_id='E1', lower_bound=E1_min, upper_bound=E1_max)
    new_E1_min, new_E1_max = 0.51, 0.52

    # act
    toy_pam.change_enzyme_bounds(enzyme_id='E1', lower_bound=new_E1_min, upper_bound=new_E1_max)
    new_toy_E1_min = toy_pam.enzyme_variables.get_by_id('E1').lower_bound
    new_toy_E1_max = toy_pam.enzyme_variables.get_by_id('E1').upper_bound

    # assert
    assert new_E1_min == new_toy_E1_min
    assert new_E1_max == new_toy_E1_max

def test_if_pamodel_change_enzyme_bounds_function_can_change_when_ub_leq_lb():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    E1_min, E1_max = 0.005, 0.5
    toy_pam.change_enzyme_bounds(enzyme_id='E1', lower_bound=E1_min, upper_bound=E1_max)
    new_E1_min, new_E1_max = 0.003, 0.004

    # act
    toy_pam.change_enzyme_bounds(enzyme_id='E1', lower_bound=new_E1_min, upper_bound=new_E1_max)
    new_toy_E1_min = toy_pam.enzyme_variables.get_by_id('E1').lower_bound
    new_toy_E1_max = toy_pam.enzyme_variables.get_by_id('E1').upper_bound

    # assert
    assert new_E1_min == new_toy_E1_min
    assert new_E1_max == new_toy_E1_max

def test_if_pamodel_sensitivity_can_be_changed_true_to_false():
    # arrange
    toy_pam = build_toy_pam(sensitivity=True)
    toy_R1_lb = -toy_pam.constraints['R1_lb'].ub
    toy_R1_ub = toy_pam.constraints['R1_ub'].ub

    # act
    toy_pam.sensitivity = False
    new_toy_R1_lb, new_toy_R1_ub = toy_pam.reactions.R1.bounds

    # assert
    assert all(true_bound == false_bound for true_bound, false_bound in zip([toy_R1_lb, toy_R1_ub],
                                                                            [new_toy_R1_lb, new_toy_R1_ub]))
    #check if constraints are removed
    assert all([constr_id not in toy_pam.constraints.keys() for constr_id in ['R1_lb', 'R1_ub']])

def test_if_pamodel_sensitivity_can_be_changed_false_to_true():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    toy_R1_lb, toy_R1_ub = toy_pam.reactions.R1.bounds

    # act
    toy_pam.sensitivity = True
    new_toy_R1_lb = -toy_pam.constraints['R1_lb'].ub
    new_toy_R1_ub = toy_pam.constraints['R1_ub'].ub

    # assert
    assert all(false_bound == true_bound for false_bound, true_bound in zip([toy_R1_lb, toy_R1_ub],
                                                                            [new_toy_R1_lb, new_toy_R1_ub]))
    #check if reaction bounds are relaxed
    assert all([toy_R1_lb>toy_pam.reactions.R1.lower_bound, toy_R1_ub<toy_pam.reactions.R1.upper_bound])


def test_if_pamodel_sensitivity_can_be_changed_true_to_false_ecolicore():
    # arrange
    ecolicore_pam = set_up_ecolicore_pam(sensitivity=True)
    glc_lb = -ecolicore_pam.constraints['EX_glc__D_e_lb'].ub
    glc_ub = ecolicore_pam.constraints['EX_glc__D_e_ub'].ub

    # act
    ecolicore_pam.sensitivity = False
    new_glc_lb, new_glc_ub = ecolicore_pam.reactions.EX_glc__D_e.bounds

    # assert
    assert all(true_bound == false_bound for true_bound, false_bound in zip([glc_lb, glc_ub],
                                                                            [new_glc_lb, new_glc_ub]))

    assert all([constr_id not in ecolicore_pam.constraints.keys() for constr_id in ['EX_glc__D_e_lb', 'EX_glc__D_e_ub']])

def test_if_pamodel_sensitivity_can_be_changed_false_to_true_ecolicore():
    # arrange
    ecolicore_pam = set_up_ecolicore_pam(sensitivity=False)
    glc_lb, glc_ub = ecolicore_pam.reactions.EX_glc__D_e.bounds

    # act
    ecolicore_pam.sensitivity = True
    new_glc_lb = -ecolicore_pam.constraints['EX_glc__D_e_lb'].ub
    new_glc_ub = ecolicore_pam.constraints['EX_glc__D_e_ub'].ub

    # assert
    assert all(false_bound == true_bound for false_bound, true_bound in zip([glc_lb, glc_ub],
                                                                            [new_glc_lb, new_glc_ub]))

def test_if_pamodel_copy_function_works():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    # act
    toy_pam_copy = toy_pam.copy()
    # assert
    assert_bounds(toy_pam, toy_pam_copy)
    assert_total_protein_content(toy_pam, toy_pam_copy)

def test_if_pamodel_copy_function_works_with_pickle():
    # arrange
    toy_pam = build_toy_pam(sensitivity=False)
    # act
    toy_pam_copy = toy_pam.copy(copy_with_pickle = True)
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

def test_if_pamodel_remove_sectors_can_remove_active_enzyme_sector():
    # Arrange
    toy_pam = build_toy_pam(sensitivity=False)

    # Act
    toy_pam.remove_sectors('ActiveEnzymeSector')

    # Assert
    assert len(toy_pam.enzymes) == 0
    assert all('EC_' not in key for key in toy_pam.constraints.keys())
    assert toy_pam.TOTAL_PROTEIN_CONSTRAINT_ID not in toy_pam.constraints.keys()

def test_if_pamodel_remove_sectors_can_remove_translational_protein_sector():
    # Arrange
    toy_pam = build_toy_pam(sensitivity=False)
    sector = toy_pam.sectors.get_by_id('TranslationalProteinSector')
    lin_rxn = toy_pam.reactions.get_by_id(sector.id_list[0])

    # get the old upperbound for the total protein constraint
    tpc_ub = toy_pam.constraints[toy_pam.TOTAL_PROTEIN_CONSTRAINT_ID].ub

    # Act
    toy_pam.remove_sectors('TranslationalProteinSector')

    #make sure the linear coefficients associated to the sector in the total protein constriant are 0
    lin_coeff = toy_pam.constraints[toy_pam.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients([
        lin_rxn.forward_variable,
        lin_rxn.reverse_variable,
    ])

    # Assert
    assert toy_pam.constraints[toy_pam.TOTAL_PROTEIN_CONSTRAINT_ID].ub == tpc_ub + sector.intercept
    assert all(coeff == 0 for coeff in lin_coeff.values())

def test_if_pamodel_gets_catalyzing_enzymes_for_enzyme_object():
    # Arrange
    sut = set_up_toy_pam_with_isozymes_and_enzymecomplex(sensitivity = False)
    enzyme_ut = 'E10'
    associated_enzymes = ['E10', 'E10_E11_E3']

    # Assert
    catalyzing_enzymes = sut._get_catalyzing_enzymes_for_enzyme(enzyme_ut)

    # Assert
    assert all(enz in catalyzing_enzymes for enz in associated_enzymes)

def test_if_rxn2protein_can_be_added_to_pam():
    # Arrange
    sut = set_up_toy_pam_with_isozymes_and_enzymecomplex(sensitivity=False)
    rxn2protein = {'R1': {'E15_E21':
                                {'f': 1,
                                'b': 0.5,
                                'molmass': 1e3,
                                'genes': ['g1', 'g2'],
                                'protein_reaction_association': [['E15_E21'], ['E11']]
                                }}}
    # Act
    sut.add_rxn2protein_to_active_enzymes(rxn2protein)

    # Arrange
    assert 'E15_E21' in sut.enzymes
    assert 'E15_E21' in sut.sectors.ActiveEnzymeSector.rxn2protein['R1'].keys()

@pytest.mark.parametrize('slope', [
    (1e-3),(-1e-3)
])
def test_if_adjust_slopes_adjust_slope_for_enzyme_sector_correctly(slope: float):
    # arrange
    sut = build_toy_pam(sensitivity=False)
    sector1 = sut.sectors.get_by_id('TranslationalProteinSector')
    sector2 = sut.sectors.get_by_id('UnusedEnzymeSector')
    sector2.slope = slope

    lin_rxn = sut.reactions.get_by_id(sector1.id_list[0])
    lin_rxn_coefficients_original = sut.constraints[sut.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients(
        [lin_rxn.forward_variable, lin_rxn.reverse_variable]
    )

    # act'
    sut._adjust_sector_slope_in_total_protein_constraint(sector2,
                                                         lin_rxn=lin_rxn)
    #make sure the linear coefficients associated to the sector in the total protein constriant are 0
    lin_coeff = sut.constraints[sut.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients([
        lin_rxn.forward_variable,
        lin_rxn.reverse_variable,
    ])


    # Assert
    assert all(coeff == pytest.approx(lin_coeff[var]-slope, rel=1e4) for var,coeff in lin_rxn_coefficients_original.items())

#######################################################################################################
#HELPER METHODS
#######################################################################################################
def build_active_enzyme_sector(Config):
    n=9
    kcat_fwd = [1, 0.5, 1, 1, 0.5 ,0.45, 1.5]  # High turnover for exchange reactions
    kcat_rev = [kcat for kcat in kcat_fwd]
    rxn2kcat = {}
    protein2gene = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id:
                                       {f'E{i + 1}':
                                            {'f': kcat_fwd[i] / (3600 * 1e-6),
                                             'b': kcat_rev[i] / (3600 * 1e-6),
                                             'molmass': 1e6,
                                             'protein_reaction_association': [[f'E{i + 1}']]}}}}
        protein2gene = {**protein2gene, **{f'E{i + 1}': [[f'g{i + 1}']]}}

    return ActiveEnzymeSector(rxn2protein=rxn2kcat, protein2gene=protein2gene, configuration=Config)

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