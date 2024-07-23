import pytest
from src.PAModelpy.configuration import Config
from src.PAModelpy.PAModel import PAModel


from Scripts.toy_ec_pam import build_toy_gem, build_active_enzyme_sector, build_translational_protein_sector, build_unused_protein_sector
from Scripts.pam_generation_uniprot_id import (set_up_ecolicore_pam, set_up_ecoli_pam, set_up_toy_pam,
                                        parse_gpr_information_for_protein2genes,
                                        parse_gpr_information_for_rxn2protein)


def test_gpr_information_is_parsed_correctly():
    # Arrange
    enzyme_id = '1.1.5.3'
    gpr_string ='(b0902 and b0903) or (b0902 and b3114) or (b3951 and b3952) or ((b0902 and b0903) and b2579)'

    parsed_gpr_reference = [['b0902','b0903'], ['b0902','b3114'], ['b3951' ,'b3952'],['b0902', 'b0903','b2579']]

    # Act
    gpr_parsed = parse_gpr_information_for_protein2genes(gpr_string)

    # Assert
    assert gpr_parsed == parsed_gpr_reference

def test_gpr_information_for_protein_is_correctly_filtered():
    # Arrange
    enzyme_id = 'P0AFH2'
    gpr_string ='((b1245 and b1247 and b1244 and b1329 and b1246) or (b2177 and b2180 and b2178 and b2179) or (b1245 and b1247 and b1246 and b1244 and b1243))'
    protein2gene = {'P02932':'b0241', 'P77348': 'b1329','P76027': 'b1246','P33913': 'b2177','P77737': 'b1247',
                    'P0AFU0': 'b2178', 'P33916': 'b2180', 'P0AFH6': 'b1245', 'P23843': 'b1243','P0AFH2': 'b1244',
                    'P33915': 'b2179'}
    gene2protein = {v: k for k, v in protein2gene.items()}

    parsed_gr_reference = [['b1245','b1247','b1244','b1329','b1246'],
                           ['b1245','b1247','b1246','b1244','b1243']]

    parsed_pr_reference = [['P0AFH6', 'P77737', 'P0AFH2', 'P77348', 'P76027'],
                           ['P0AFH6', 'P77737', 'P76027', 'P0AFH2', 'P23843']]

    # Act
    gr_parsed, pr_parsed = parse_gpr_information_for_rxn2protein(gpr_string, gene2protein, protein2gene, enzyme_id)

    # Assert
    assert gr_parsed == parsed_gr_reference
    assert pr_parsed == parsed_pr_reference

def test_if_enzyme_complex_in_toy_pam_is_parsed_correctly():
    sut = set_up_toy_pam_with_enzyme_complex(sensitivity=False)

    assert all([enz in sut.enzymes for enz in ['E1', 'E2', 'E10', 'E2_E10']])
    assert all([const not in sut.constraints.keys() for const in ['EC_E10_f', 'EC_E2_f']])
    constraint = sut.constraints['EC_E2_E10_f'].get_linear_coefficients([sut.reactions.CE_R2_E2_E10.forward_variable])
    assert constraint[sut.reactions.CE_R2_E2_E10.forward_variable] > 0

def test_if_isozymes_in_toy_pam_are_parsed_correctly():
    sut = set_up_toy_pam_with_isozymes(sensitivity=False)

    #check whether catalytic reactions are configured correclty
    neg_constraint = sut.constraints['CE_R2'].get_linear_coefficients([sut.reactions.R2.reverse_variable,
                                                                   sut.reactions.CE_R2_E2.forward_variable,
                                                                   sut.reactions.CE_R2_E10.forward_variable])

    pos_constraint = sut.constraints['CE_R2'].get_linear_coefficients([sut.reactions.R2.forward_variable,
                                                                   sut.reactions.CE_R2_E2.reverse_variable,
                                                                   sut.reactions.CE_R2_E10.reverse_variable])

    assert all([enz in sut.enzymes for enz in ['E1', 'E2', 'E10']])
    assert all([const in sut.constraints.keys() for const in ['EC_E10_f', 'EC_E2_f']])
    assert all([rxn in sut.reactions for rxn in ['CE_R1_E1', 'CE_R2_E2', 'CE_R2_E10']])
    assert all([v == -1 for v in neg_constraint.values()])
    assert all([v == 1 for v in pos_constraint.values()])

def test_if_toy_pam_with_isozymes_has_same_growth_rate_as_without():
    sut = set_up_toy_pam_with_isozymes(sensitivity=False)
    toy_pam = set_up_toy_pam()

    for model in [sut, toy_pam]:
        model.change_reaction_bounds('R1', 0,0.01)
        model.optimize()

    assert sut.objective.value == pytest.approx(toy_pam.objective.value, abs = 1e-6)

def test_if_toy_pam_with_isozymes_and_enzymecomplex_has_same_growth_rate_as_without():
    sut = set_up_toy_pam_with_isozymes_and_enzymecomplex(sensitivity=False)
    toy_pam = set_up_toy_pam()

    for model in [sut, toy_pam]:
        model.change_reaction_bounds('R1', 0,0.01)
        model.optimize()

    assert sut.objective.value == pytest.approx(toy_pam.objective.value, abs = 1e-6)


def test_if_toy_pam_with_enzyme_comples_has_same_growth_rate_as_without():
    sut = set_up_toy_pam_with_enzyme_complex(sensitivity=False)
    toy_pam = set_up_toy_pam()

    for model in [sut, toy_pam]:
        model.change_reaction_bounds('R1', 0,0.01)
        model.optimize()

    assert sut.objective.value == pytest.approx(toy_pam.objective.value, abs = 1e-6)

# def test_set_up_ecolicore_pam_works():
#     sut = set_up_ecolicore_pam()
#     assert True
# def test_set_up_ecoli_pam_works():
#     sut = set_up_ecoli_pam()
#     assert True

#########################################################################################################################
# HELPER FUNCTIONS
##########################################################################################################################

def set_up_toy_pam_with_enzyme_complex(sensitivity =True):
    config = Config()
    #setting the configuration for the toy model
    config.BIOMASS_REACTION = 'R7'
    config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    config.CO2_EXHANGE_RXNID = 'R8'
    config.ACETATE_EXCRETION_RXNID = 'R9'

    Etot = 0.6*1e-3
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(config)
    #add an enzyme associated to enzyme complex to the toy model
    active_enzyme.rxn2protein['R2']['E2']['protein_reaction_association'] = [['E2', 'E10']]
    active_enzyme.rxn2protein['R2']['E10']= active_enzyme.rxn2protein['R2']['E2'].copy()

    #build the toy model
    unused_enzyme = build_unused_protein_sector(config)
    translation_enzyme = build_translational_protein_sector(config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme, p_tot=Etot, sensitivity=sensitivity)
    pamodel.objective = 'R7'
    config.reset()
    return pamodel

def set_up_toy_pam_with_isozymes(sensitivity =True):
    config = Config()
    #setting the configuration for the toy model
    config.BIOMASS_REACTION = 'R7'
    config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    config.CO2_EXHANGE_RXNID = 'R8'
    config.ACETATE_EXCRETION_RXNID = 'R9'

    Etot = 0.6*1e-3
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(config)
    #add an enzyme associated to isozymes to the toy model
    active_enzyme.rxn2protein['R2']['E2']['protein_reaction_association'] = [['E2'], ['E10']]
    active_enzyme.rxn2protein['R2']['E10']= active_enzyme.rxn2protein['R2']['E2'].copy()

    #build the toy model
    unused_enzyme = build_unused_protein_sector(config)
    translation_enzyme = build_translational_protein_sector(config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme, p_tot=Etot, sensitivity=sensitivity)
    pamodel.objective = 'R7'
    config.reset()
    return pamodel

def set_up_toy_pam_with_isozymes_and_enzymecomplex(sensitivity =True):
    config = Config()
    #setting the configuration for the toy model
    config.BIOMASS_REACTION = 'R7'
    config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    config.CO2_EXHANGE_RXNID = 'R8'
    config.ACETATE_EXCRETION_RXNID = 'R9'

    Etot = 0.6*1e-3
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(config)
    #add an enzyme associated to isozymes to the toy model
    active_enzyme.rxn2protein['R2']['E2']['protein_reaction_association'] = [['E2'], ['E10']]
    active_enzyme.rxn2protein['R2']['E10']= active_enzyme.rxn2protein['R2']['E2'].copy()

    #add an enzyme associated to enzyme complex to the toy model
    active_enzyme.rxn2protein['R3']['E3']['protein_reaction_association'] = [['E3','E10', 'E11']]
    active_enzyme.rxn2protein['R3']['E10']= active_enzyme.rxn2protein['R3']['E3'].copy()
    active_enzyme.rxn2protein['R3']['E11']= active_enzyme.rxn2protein['R3']['E3'].copy()


    #build the toy model
    unused_enzyme = build_unused_protein_sector(config)
    translation_enzyme = build_translational_protein_sector(config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme, p_tot=Etot, sensitivity=sensitivity)
    pamodel.objective = 'R7'
    config.reset()
    return pamodel
