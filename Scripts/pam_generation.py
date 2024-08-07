import cobra
import pandas as pd
import os
from typing import Union

# load PAMpy modules
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.PAModelpy.configuration import Config

from src.PAModelpy import EnzymeVariable

from Scripts.toy_ec_pam import build_toy_gem, build_active_enzyme_sector, build_translational_protein_sector, build_unused_protein_sector

'Function library for making Protein Allocation Models as described in the publication'


def set_up_toy_pam(sensitivity =True):
    config = Config()
    #setting the configuration for the toy model
    config.BIOMASS_REACTION = 'R7'
    config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    config.CO2_EXHANGE_RXNID = 'R8'
    config.ACETATE_EXCRETION_RXNID = 'R9'

    Etot = 0.6*1e-3
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(config)
    unused_enzyme = build_unused_protein_sector(config)
    translation_enzyme = build_translational_protein_sector(config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme, p_tot=Etot, sensitivity=sensitivity)
    pamodel.objective = 'R7'
    config.reset()
    return pamodel

def set_up_ecolicore_pam(total_protein:bool = True, active_enzymes: bool = True, translational_enzymes:bool = True, unused_enzymes:bool = True, sensitivity =True):
    # Setting the relative paths
    PAM_DATA_FILE_PATH = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')

    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]

    # load the genome-scale information
    model = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))

    #load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        enzyme_db = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes')
        enzyme_db = enzyme_db.set_index('rxnID')
        # correct reaction IDS
        for idx in enzyme_db.index.to_list():
            # transprt reactions<

            if 'pp' in idx:
                idx_new = idx.replace('pp', '')
                if idx_new not in enzyme_db.index:
                    enzyme_db.rename(index={idx: idx_new}, inplace=True)
            if 'ex' in idx:
                idx_new = idx.replace('ex', '')
                if idx_new not in enzyme_db.index:
                    enzyme_db.rename(index={idx: idx_new}, inplace=True)

                    # replace NaN values with unique identifiers
        # replace NaN enzyme ids with a dummy enzyme identifier
        # select the NaN values
        nan_values = enzyme_db['EC_nmbr'].isnull()
        # make a list with unique ids
        nan_ids = [f'E{i}' for i in range(nan_values.sum())]
        # replace nan values by unique id
        enzyme_db.loc[nan_values, 'EC_nmbr'] = nan_ids

        # create enzyme objects for each gene-associated reaction
        kcats = {}
        rxn2ec = {}
        molmass = {}
        for rxn in model.reactions:
            if rxn.genes:
                # correct transport reactions
                if 't' in rxn.id:
                    rxn.id = rxn.id
                # are enzyme information in the PAM database?
                rev = 0  # denotes reversibility
                if rxn.lower_bound >= 0:
                    # irreversible reaction (forward direction)
                    rev = 0
                    rxn_id = rxn.id  # save reaction ID for retrieveing molar masses/enzyme information later
                    if rxn.id in enzyme_db.index:
                        kcats[rxn.id] = {'f': enzyme_db.loc[rxn.id, 'kcat']}
                elif rxn.upper_bound <= 0:
                    # irreversible reaction (reverse direction)
                    rev = 1
                    rxn_id = rxn.id + '_b'
                    if rxn_id in enzyme_db.index:
                        kcats[rxn.id] = {'b': enzyme_db.loc[rxn_id, 'kcat']}
                else:
                    rev = 2
                    # reversible reaction
                    rxn_id_f = rxn.id + '_f'
                    rxn_id_b = rxn.id + '_b'
                    if rxn_id_f in enzyme_db.index and rxn_id_b in enzyme_db.index:
                        rxn_id = rxn_id_f  # save reaction ID for retrieveing molar masses/enzyme information later
                        kcats[rxn.id] = {'f': enzyme_db.loc[rxn_id_f, 'kcat'],
                                         'b': enzyme_db.loc[rxn_id_b, 'kcat']}

                    else:
                        # try if only forward reaction is in database
                        rxn_id = rxn.id  # save reaction ID for retrieveing molar masses/enzyme information later
                        kcats[rxn.id] = {'f': enzyme_db.loc[rxn.id, 'kcat'],
                                         'b': enzyme_db.loc[
                                                  rxn.id, 'kcat'] / 2}  # deduce backwards kcat from forward value

                # where enzyme information found?
                if rxn.id in kcats.keys():
                    # save molmass
                    molmass[rxn.id] = enzyme_db.loc[rxn_id, 'molMass']
                    # save enzyme information
                    # is enzyme information NaN?
                    if pd.isna(enzyme_db.loc[rxn_id, 'EC_nmbr']):
                        rxn2ec[rxn.id] = ''
                    else:
                        rxn2ec[rxn.id] = enzyme_db.loc[rxn_id, 'EC_nmbr']


                else:
                    # no enzyme information found
                    print('No enzyme information found for reaction: ' + rxn.id)
                    # Create generic Enzyme with mean molar masses and kcat
                    if rev == 0:
                        kcats[rxn.id] = {'f': 22}
                    elif rev == 1:
                        kcats[rxn.id] = {'b': 22}
                    else:
                        kcats[rxn.id] = {'f': 22, 'b': 22}

                    molmass[rxn.id] = 3.947778784340140e04

        rxn2protein = {}
        for rxn, ec in rxn2ec.items():
            ec_dict = {**kcats[rxn], **{'molmass': molmass[rxn]}}
            # add enzyme to enzymes related to reaction if these are already stored
            if rxn in rxn2protein.keys():
                rxn2protein[rxn] = {**rxn2protein[rxn], **{ec: ec_dict}}
            # if not create new reaction entry
            else:
                rxn2protein[rxn] = {ec: ec_dict}

        # create active enzymes sector
        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        # translational protein sector parameter (substrate dependent)
        id_list_tps = ['EX_glc__D_e']
        tps_0 = [0.04992]  # g/gDW
        tps_mu = [-0.002944]  # g h/gDW -> transformed to match glucose uptake variable
        molmass_tps = [405903.94]  # g/mol

        # translational protein sector
        translation_enzyme_sector = TransEnzymeSector(
            id_list=id_list_tps,
            tps_0=tps_0,
            tps_mu=tps_mu,
            mol_mass=molmass_tps,
        )
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        id_list_ups = [BIOMASS_REACTION]
        ups_0 = [0.0407]  # g/gDW
        ups_mu = [-0.0214]  # g h/gDW -> negative relation with growth rate
        molmass_ups = [405903.94]  # g/mol

        unused_enzyme_sector = UnusedEnzymeSector(
            id_list=id_list_ups,
            ups_0=ups_0,
            ups_mu=ups_mu,
            mol_mass=molmass_ups,
        )
    else:
        unused_enzyme_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pa_model = PAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector, unused_sector=unused_enzyme_sector)
    return pa_model

def set_up_ecoli_pam(total_protein: Union[bool, float] = True, active_enzymes: bool = True,
                   translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True):

    config = Config()
    config.reset()
    # Setting the relative paths
    pam_info_file = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')

    # some other constants
    TOTAL_PROTEIN_CONCENTRATION = 0.258  # [g_prot/g_cdw]

    #setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join('Models', 'iML1515.xml'))

    #check if a different total protein concentration is given
    if isinstance(total_protein, float):
        TOTAL_PROTEIN_CONCENTRATION = total_protein

    # load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        active_enzyme_info_old = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')

        # replace NaN values with unique identifiers
        # select the NaN values
        nan_values = active_enzyme_info_old['EC_nmbr'].isnull()
        # make a list with unique ids
        nan_ids = [f'E{i}' for i in range(nan_values.sum())]
        # replace nan values by unique id
        active_enzyme_info_old.loc[nan_values, 'EC_nmbr'] = nan_ids

        # parse the enzyme information (kcat values, identifiers and molmasses)
        kcats_dict = active_enzyme_info_old.set_index(keys='rxnID').loc[:, 'kcat'].to_dict()
        ec_dict = active_enzyme_info_old.set_index(keys='rxnID').loc[:, 'EC_nmbr'].to_dict()
        molmass_dict = mol_mass = active_enzyme_info_old.set_index(keys='rxnID').loc[:, 'molMass'].to_dict()

        kcats = {}
        # save fwd and bckw kcats separately in the form of: {rxn_id: {'f': kcat_f, 'b': kcat_b}}
        for rxn, kcat in kcats_dict.items():
            # reversible reaction
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                direction = rxn[-1]
                # check if the reaction already exists in the kcat dictionary
                try:
                    kcats[rxn[:-2]][direction] = kcat
                except:
                    kcats[rxn[:-2]] = {direction: kcat}
            # irreversible reaction
            else:
                kcats[rxn] = {'f': kcat}

        rxn2ec = {}
        # parse the enzyme identifiers for the reactions
        for rxn, ec in ec_dict.items():
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                rxn = rxn[:-2]
            for enz in str(ec).split(','):
                rxn2ec[rxn] = enz.strip()

        molmass = {}
        # parse the enzyme molmasses for the reactions
        for rxn, mw in molmass_dict.items():
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                rxn = rxn[:-2]
            molmass[rxn] = mw

        rxn2protein = {}
        for rxn, ec in rxn2ec.items():
            ec_dict = {**kcats[rxn], **{'molmass': molmass[rxn]}}
            # add enzyme to enzymes related to reaction if these are already stored
            if rxn in rxn2protein.keys():
                rxn2protein[rxn] = {**rxn2protein[rxn], **{ec: ec_dict}}
            # if not create new reaction entry
            else:
                rxn2protein[rxn] = {ec: ec_dict}

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein, configuration = config)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        translational_info = pd.read_excel(pam_info_file, sheet_name='Translational')
        translation_enzyme_sector = TransEnzymeSector(
            id_list=[translational_info[translational_info.Parameter == 'id_list'].loc[0, 'Value']],
            tps_0=[translational_info[translational_info.Parameter == 'tps_0'].loc[1, 'Value']],
            tps_mu=[-translational_info[translational_info.Parameter == 'tps_mu'].loc[2, 'Value']],
            mol_mass=[translational_info[translational_info.Parameter == 'mol_mass'].loc[3, 'Value']],
            configuration = config)
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        unused_protein_info = pd.read_excel(pam_info_file, sheet_name='ExcessEnzymes')

        ups_0 = unused_protein_info[unused_protein_info.Parameter == 'ups_0'].loc[2, 'Value']
        smax = unused_protein_info[unused_protein_info.Parameter == 's_max_uptake'].loc[1, 'Value']

        unused_protein_sector = UnusedEnzymeSector(
            id_list=[unused_protein_info[unused_protein_info.Parameter == 'id_list'].loc[0, 'Value']],
            ups_mu=[ups_0 / smax],
            ups_0=[ups_0],
            mol_mass=[unused_protein_info[unused_protein_info.Parameter == 'mol_mass'].loc[3, 'Value']],
            configuration = config)
    else:
        unused_protein_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pamodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_protein_sector, sensitivity=sensitivity, configuration = config
                      )
    return pamodel

def parse_coefficients(pamodel):
    Ccsc = list()

    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
        Ccsc += pamodel.capacity_sensitivity_coefficients[
            pamodel.capacity_sensitivity_coefficients['constraint'] == csc].coefficient.to_list()

    Cesc = pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

    return Ccsc, Cesc

def parse_esc(pamodel):
    return pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

if __name__ == '__main__':
    ecoli_pam = set_up_ecoli_pam(sensitivity=False)
    # ecoli_pam.objective = ecoli_pam.BIOMASS_REACTION
    ecoli_pam.change_reaction_bounds('EX_glc__D_e', -10, 0)
    ecoli_pam.optimize()
    print(ecoli_pam.objective.value)
    import pickle

    with open('path_to_your_pickle_file.pkl', 'wb') as file:
        p = pickle.dump(ecoli_pam, file)
    with open('path_to_your_pickle_file.pkl', 'rb') as file:
        ob = pickle.load(file)