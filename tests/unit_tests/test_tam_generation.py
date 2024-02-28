import pytest
import os
import pandas as pd
import cobra

from src.PAModelpy import PAModel, UnusedEnzymeSector, ActiveEnzymeSector, TransEnzymeSector
from Scripts.tam_generation import set_up_toy_tam, set_up_ecolicore_tam
from Scripts.pam_generation import set_up_ecolicore_pam, set_up_toy_pam


# def test_set_up_ecolicore_pam_builds_correct_model():
#     # Arrange
#     sut = set_up_ecolicore_pam()
#     ecolicore_pam_original = setup_ecolicore_pam_validation()
#
#     print(ecolicore_pam_original.enzymes)
#     # print(sut.enzymes)
#     # Assert
#     for enzyme in sut.enzyme_variables:
#         if enzyme.id not in ecolicore_pam_original.enzymes:
#             print(enzyme.id, enzyme.reactions)
#
#     for rxn in sut.reactions:
#         if rxn.id not in ecolicore_pam_original.reactions:
#             print(rxn.id)
#
#     assert False

def test_set_up_toy_tam_works():
    sut = set_up_toy_tam

def test_set_up_ecolicore_tam_works():
    sut = set_up_ecolicore_tam()
    assert True

def test_toy_tam_without_constraints_has_same_result_as_pam():
    # Arrange
    sut = set_up_toy_tam()
    toy_pam = set_up_toy_pam()

    sut.change_reaction_bounds(rxn_id='R1',
                                   lower_bound=0, upper_bound=1e-3)
    toy_pam.change_reaction_bounds(rxn_id='R1',
                                   lower_bound=0, upper_bound=1e-3)


    # Act
    sut.optimize()
    toy_pam.optimize()

    # Arrange
    assert 'optimal' == toy_pam.solver.status
    assert 'optimal' == sut.solver.status
    assert toy_pam.objective.value == sut.objective.value

def test_ecolicore_tam_without_constraints_has_same_result_as_pam():
    # Arrange
    sut = set_up_ecolicore_tam()
    ecolicore_pam = set_up_ecolicore_pam()

    sut.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                   lower_bound=-10, upper_bound=0)
    ecolicore_pam.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                   lower_bound=-10, upper_bound=0)

    # Act
    sut.optimize()
    ecolicore_pam.optimize()
    # for enzyme in ecolicore_pam.enzyme_variables:
    #     print(enzyme.id, enzyme.flux)

    # Arrange
    assert 'optimal' == ecolicore_pam.solver.status
    assert 'optimal' == sut.solver.status
    assert_solution_pam_equals_tam(ecolicore_pam, sut)
    assert ecolicore_pam.objective.value == sut.objective.value

def assert_solution_pam_equals_tam(pam, tam):
    # for reaction in tam.reactions:
    #     # print(reaction.id, pam.reactions.get_by_id(reaction.id).flux,reaction.flux)
    #     # assert pam.reactions.get_by_id(reaction.id).flux == pytest.approx(reaction.flux, rel= 1e-3)
    # for enzyme in tam.enzyme_variables:
    #     if enzyme.id == 'E114':
            # print(enzyme.id,enzyme.reactions, pam.enzyme_variables.get_by_id(enzyme.id).flux,enzyme.flux)

        # assert pam.enzymes.get_by_id(enzyme.id).flux == pytest.approx(enzyme.flux, rel= 1e-3)
    assert pam.objective.value == tam.objective.value



# def setup_ecolicore_pam_validation(total_protein:bool = True, active_enzymes: bool = True, translational_enzymes:bool = True, unused_enzymes:bool = True, sensitivity =True):
#     # Setting the relative paths
#     DATA_DIR = os.path.join(os.getcwd(), 'Data')
#     MODEL_DIR = os.path.join(os.getcwd(), 'Models')
#     PAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')
#
#     # PAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'proteinAllocationModel_iML1515_EnzymaticData_230503.xls')
#
#     # some other constants
#     BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
#     TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]
#
#     # load the genome-scale information
#     model = cobra.io.load_json_model(os.path.join(MODEL_DIR, 'e_coli_core.json'))
#
#     #load example data for the E.coli iML1515 model
#     if active_enzymes:
#         # load active enzyme sector information
#         enzyme_db = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes')
#         enzyme_db = enzyme_db.set_index('rxnID')
#         # correct reaction IDS
#         for idx in enzyme_db.index.to_list():
#             # transprt reactions<
#
#             if 'pp' in idx:
#                 idx_new = idx.replace('pp', '')
#                 if idx_new not in enzyme_db.index:
#                     enzyme_db.rename(index={idx: idx_new}, inplace=True)
#             if 'ex' in idx:
#                 idx_new = idx.replace('ex', '')
#                 if idx_new not in enzyme_db.index:
#                     enzyme_db.rename(index={idx: idx_new}, inplace=True)
#
#                     # replace NaN values with unique identifiers
#         # replace NaN enzyme ids with a dummy enzyme identifier
#         # select the NaN values
#         nan_values = enzyme_db['EC_nmbr'].isnull()
#         # make a list with unique ids
#         nan_ids = [f'E{i}' for i in range(nan_values.sum())]
#         # replace nan values by unique id
#         enzyme_db.loc[nan_values, 'EC_nmbr'] = nan_ids
#
#         # create enzyme objects for each gene-associated reaction
#
#         # parse the enzyme information (kcat values, identifiers and molmasses)
#         # create enzyme objects for each gene-associated reaction
#         kcats = {}
#         rxn2ec = {}
#         molmass = {}
#         for rxn in model.reactions:
#             if rxn.genes:
#                 # correct transport reactions
#                 if 't' in rxn.id:
#                     rxn.id = rxn.id
#                 # are enzyme information in the PAM database?
#                 rev = 0  # denotes reversibility
#                 if rxn.lower_bound >= 0:
#                     # irreversible reaction (forward direction)
#                     rev = 0
#                     rxn_id = rxn.id  # save reaction ID for retrieveing molar masses/enzyme information later
#                     if rxn.id in enzyme_db.index:
#                         kcats[rxn.id] = {'f': enzyme_db.loc[rxn.id, 'kcat']}
#                 elif rxn.upper_bound <= 0:
#                     # irreversible reaction (reverse direction)
#                     rev = 1
#                     rxn_id = rxn.id + '_b'
#                     if rxn_id in enzyme_db.index:
#                         kcats[rxn.id] = {'b': enzyme_db.loc[rxn_id, 'kcat']}
#                 else:
#                     rev = 2
#                     # reversible reaction
#                     rxn_id_f = rxn.id + '_f'
#                     rxn_id_b = rxn.id + '_b'
#                     if rxn_id_f in enzyme_db.index and rxn_id_b in enzyme_db.index:
#                         rxn_id = rxn_id_f  # save reaction ID for retrieveing molar masses/enzyme information later
#                         kcats[rxn.id] = {'f': enzyme_db.loc[rxn_id_f, 'kcat'],
#                                          'b': enzyme_db.loc[rxn_id_b, 'kcat']}
#
#                     else:
#                         # try if only forward reaction is in database
#                         rxn_id = rxn.id  # save reaction ID for retrieveing molar masses/enzyme information later
#                         kcats[rxn.id] = {'f': enzyme_db.loc[rxn.id, 'kcat'],
#                                          'b': enzyme_db.loc[
#                                                   rxn.id, 'kcat'] / 2}  # deduce backwards kcat from forward value
#
#                 # where enzyme information found?
#                 if rxn.id in kcats.keys():
#                     # save molmass
#                     molmass[rxn.id] = enzyme_db.loc[rxn_id, 'molMass']
#                     # save enzyme information
#                     # is enzyme information NaN?
#                     if pd.isna(enzyme_db.loc[rxn_id, 'EC_nmbr']):
#                         rxn2ec[rxn.id] = ''
#                     else:
#                         rxn2ec[rxn.id] = enzyme_db.loc[rxn_id, 'EC_nmbr']
#
#
#                 else:
#                     # no enzyme information found
#                     print('No enzyme information found for reaction: ' + rxn.id)
#                     # Create generic Enzyme with mean molar masses and kcat
#                     if rev == 0:
#                         kcats[rxn.id] = {'f': 22}
#                     elif rev == 1:
#                         kcats[rxn.id] = {'b': 22}
#                     else:
#                         kcats[rxn.id] = {'f': 22, 'b': 22}
#
#                     molmass[rxn.id] = 3.947778784340140e04
#
#         rxn2protein = {}
#         for rxn, ec in rxn2ec.items():
#             ec_dict = {**kcats[rxn], **{'molmass': molmass[rxn]}}
#             # add enzyme to enzymes related to reaction if these are already stored
#             if rxn in rxn2protein.keys():
#                 rxn2protein[rxn] = {**rxn2protein[rxn], **{ec: ec_dict}}
#             # if not create new reaction entry
#             else:
#                 rxn2protein[rxn] = {ec: ec_dict}
#
#         # create active enzymes sector
#         active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein)
#
#     else:
#         active_enzyme_sector = None
#
#     if translational_enzymes:
#         # translational protein sector parameter (substrate dependent)
#         id_list_tps = ['EX_glc__D_e']
#         tps_0 = [0.04992]  # g/gDW
#         tps_mu = [-0.002944]  # g h/gDW -> transformed to match glucose uptake variable
#         molmass_tps = [405903.94]  # g/mol
#
#         # translational protein sector
#         translation_enzyme_sector = TransEnzymeSector(
#             id_list=id_list_tps,
#             tps_0=tps_0,
#             tps_mu=tps_mu,
#             mol_mass=molmass_tps,
#         )
#     else:
#         translation_enzyme_sector = None
#
#     if unused_enzymes:
#         id_list_ups = [BIOMASS_REACTION]
#         ups_0 = [0.0407]  # g/gDW
#         ups_mu = [-0.0214]  # g h/gDW -> negative relation with growth rate
#         molmass_ups = [405903.94]  # g/mol
#
#         unused_enzyme_sector = UnusedEnzymeSector(
#             id_list=id_list_ups,
#             ups_0=ups_0,
#             ups_mu=ups_mu,
#             mol_mass=molmass_ups,
#         )
#     else:
#         unused_enzyme_sector = None
#
#     if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION
#
#     pa_model = PAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity,
#                        active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector, unused_sector=unused_enzyme_sector)
#     return pa_model