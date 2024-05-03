import cobra
import pandas as pd
import os
from typing import Union
import ast
import numpy as np

# load PAMpy modules
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.PAModelpy.configuration import Config

from src.TAModelpy import TAModel, ActivemRNASector
from Scripts.toy_ec_tam import build_toy_gem, build_active_enzyme_sector, build_translational_protein_sector, build_unused_protein_sector

'Function library for making Protein Allocation Models as described in the publication'


def set_up_toy_tam(sensitivity =True):
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

def set_up_ecolicore_tam(total_protein:bool = True, total_mrna:bool =True,active_enzymes: bool = True,
                         translational_enzymes:bool = True, unused_enzymes:bool = True,
                         sensitivity =True):
    # Setting the relative paths
    DATA_DIR = os.path.join('Data')
    MODEL_DIR = os.path.join('Models')
    TAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'TAModel','2024-04-24_gene_enzyme_reaction_relation_Ecoli.xlsx')


    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    # TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]
    TOTAL_PROTEIN_CONCENTRATION = 0.185  # [g_prot/g_cdw]

    MRNA_MU = 0.00013049558330984208 # [g_mrna/g_cdw/h]
    MRNA_0= 1.7750480089801658e-05 # [g_mrna/g_cdw]

    # MRNA_MU = 0.0036378316565569466
    # MRNA_0= 10000000
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = BIOMASS_REACTION

    # load the genome-scale information
    model = cobra.io.load_json_model(os.path.join(MODEL_DIR, 'e_coli_core.json'))

    #load example data for the E.coli iML1515 model
    if active_enzymes:
        enzyme_db = _parse_enzyme_information_from_file(TAM_DATA_FILE_PATH)

        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene, gene2transcript = parse_reaction2protein2gene2transcript(enzyme_db, model)

        # import random
        # gene2transcript2 ={}
        # for i in range(50):
        #     gene, transcript = random.choice(list(gene2transcript.items()))
        #     gene2transcript2 = {**gene2transcript2, gene:transcript}
        # gene2, transcript2 = random.choice(list(gene2transcript.items()))
        # gene3, transcript3 = random.choice(list(gene2transcript.items()))
        #
        # gene2transcript = {gene:transcript, gene2:transcript2, gene3:transcript3}
        # create active enzymes sector
        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein,
                                                  protein2gene = protein2gene,
                                                  configuration = config)

        # gene2transcript = {'gene_dummy': {'id': 'mRNA_gene_dummy', 'length': 750.0}}
        active_mrna_sector = ActivemRNASector(mrnas_0 = MRNA_0,
                                              mrnas_mu = MRNA_MU,
                                              id_list = [BIOMASS_REACTION],
                                              gene2transcript = gene2transcript,
                                              configuration = config)

    else:
        active_enzyme_sector = None
        active_mrna_sector = None

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
            configuration = config
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
            configuration = config
        )
    else:
        unused_enzyme_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    ta_model = TAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity, mrna_sector = active_mrna_sector,
            configuration = config, active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_enzyme_sector, total_mrna_constraint = total_mrna)
    return ta_model

def set_up_ecoli_tam(total_protein: Union[bool, float] = True, total_mrna:bool =True, active_enzymes: bool = True,
                   translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True):

    config = Config()
    config.reset()
    # Setting the relative paths
    MODEL_DIR = 'Models'
    TAM_DATA_FILE_PATH = os.path.join('Data', 'TAModel','2024-04-24_gene_enzyme_reaction_relation_Ecoli.xlsx')
    pam_info_file = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')

    # some other constants
    TOTAL_PROTEIN_CONCENTRATION = 0.258  # [g_prot/g_cdw]

    MRNA_MU = 0.0019274475344656643 # [g_mrna/g_cdw/h]
    MRNA_0= 0.00026217836816312917 # [g_mrna/g_cdw]

    #setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join(MODEL_DIR, 'iML1515.xml'))

    #check if a different total protein concentration is given
    if isinstance(total_protein, float):
        TOTAL_PROTEIN_CONCENTRATION = total_protein

    config = Config()
    config.reset()

    # load example data for the E.coli iML1515 model
    if active_enzymes:
        enzyme_db = _parse_enzyme_information_from_file(TAM_DATA_FILE_PATH)

        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene, gene2transcript = parse_reaction2protein2gene2transcript(enzyme_db, model)
        # create active enzymes sector
        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein,
                                                  protein2gene=protein2gene,
                                                  configuration=config)

        # gene2transcript = {'gene_dummy': {'id': 'mRNA_gene_dummy', 'length': 750.0}}
        active_mrna_sector = ActivemRNASector(mrnas_0=MRNA_0,
                                              mrnas_mu=MRNA_MU,
                                              id_list=[config.BIOMASS_REACTION],
                                              gene2transcript=gene2transcript,
                                              configuration=config)

    else:
        active_enzyme_sector = None
        active_mrna_sector = None

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

        unused_enzyme_sector = UnusedEnzymeSector(
            id_list=[unused_protein_info[unused_protein_info.Parameter == 'id_list'].loc[0, 'Value']],
            ups_mu=[ups_0 / smax],
            ups_0=[ups_0],
            mol_mass=[unused_protein_info[unused_protein_info.Parameter == 'mol_mass'].loc[3, 'Value']],
            configuration = config)
    else:
        unused_enzyme_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    ta_model = TAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity, mrna_sector=active_mrna_sector,
                       configuration=config, active_sector=active_enzyme_sector,
                       translational_sector=translation_enzyme_sector,
                       unused_sector=unused_enzyme_sector, total_mrna_constraint=total_mrna)
    return ta_model

def parse_coefficients(pamodel):
    Ccsc = list()

    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
        Ccsc += pamodel.capacity_sensitivity_coefficients[
            pamodel.capacity_sensitivity_coefficients['constraint'] == csc].coefficient.to_list()

    Cesc = pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

    return Ccsc, Cesc

def parse_esc(pamodel):
    return pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

def _parse_enzyme_information_from_file(file_path:str):
    # load active enzyme sector information
    enzyme_db = pd.read_excel(file_path, sheet_name='enzyme-gene-reaction')
    # enzyme_db = enzyme_db.set_index('rxnID')
    # correct reaction IDS
    for idx in enzyme_db.rxnID.to_list():
        # transport reactions
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

    return enzyme_db


def parse_reaction2protein2gene2transcript(enzyme_db: pd.DataFrame, model:cobra.Model) -> dict:
    # Initialize dictionaries
    rxn2protein = {}
    protein2gene = {}
    gene2transcript = {'gene_dummy': {'id': 'mRNA_gene_dummy', 'length': 750.0}}

    # Iterate over each row in the DataFrame
    for index, row in enzyme_db.iterrows():
        # Parse data from the row
        rxn_id = row['rxnID']
        rxn_id = _check_rxn_identifier_format(rxn_id)
        #only parse those reactions which are in the model
        kcat_f_b = _get_fwd_bckw_kcat(rxn_id, row['kcat'], model)
        if kcat_f_b is None: continue
        #check if there is a forward or backward string which needs to be removed from the rxn id
        if any([rxn_id.endswith('_f'), rxn_id.endswith('_b')]): rxn_id = rxn_id[:-2]
        rxn = model.reactions.get_by_id(rxn_id)
        #get the identifiers and replace nan values by dummy placeholders
        enzyme_id = row['EC_nmbr']
        gene_id = row['gene_id']
        if len(rxn.genes)>0:
            if not isinstance(enzyme_id, str): enzyme_id = 'Enzyme_'+rxn_id
            if not isinstance(gene_id, str): gene_id = 'gene_dummy'
            mrna_length = row['mrna_length']

            #get the gene-protein-reaction-associations
            gpr = parse_gpr_information_for_protein2genes(row['gpr'], enzyme_id)
            # Create gene-protein-reaction associations
            if enzyme_id not in protein2gene:
                protein2gene[enzyme_id] = gpr
            if enzyme_id in protein2gene:
                #assume that there is a single protein if the previous relation was not assigned a gene
                if 'gene_' in protein2gene[enzyme_id][0][0]:
                    protein2gene[enzyme_id] = gpr
                elif 'gene_' in gpr[0][0]:
                    continue
                #check if there are multiple enzymes for a single EC number by checking if the gpr entries are different
                #get all entries with this enzyme_id:
                protein2gene_entries = [key for key in protein2gene.keys() if enzyme_id in key]

                if not any([sorted([sorted(rel) for rel in protein2gene[protein]]) == sorted([sorted(rel) for rel in gpr]) for protein in protein2gene_entries]):
                    protein2gene[f'{enzyme_id}_{rxn_id}'] = gpr
                    if rxn_id not in rxn2protein.keys():
                        rxn2protein[rxn_id]= {f'{enzyme_id}_{rxn_id}':{
                            'f': kcat_f_b[0],  # Forward kcat
                            'b': kcat_f_b[1],  # Backward kcat
                            'molmass': row['molMass'],
                            'genes': gpr,
                            'complex_with': None
                        }}
                    else:
                        rxn2protein[rxn_id][f'{enzyme_id}_{rxn_id}']= {
                            'f': kcat_f_b[0],  # Forward kcat
                            'b': kcat_f_b[1],  # Backward kcat
                            'molmass': row['molMass'],
                            'genes': gpr,
                            'complex_with': None
                        }

            # Create gene2transcript dictionary
            if gene_id not in gene2transcript.keys():
                gene2transcript[gene_id] = {'id': f'mRNA_{gene_id}', 'length': mrna_length}

        # Create rxn2protein dictionary
        if rxn_id not in rxn2protein:
            rxn2protein[rxn_id] = {}
        if enzyme_id not in rxn2protein[rxn_id]:
            rxn2protein[rxn_id][enzyme_id] = {
                'f': kcat_f_b[0],  # Forward kcat
                'b': kcat_f_b[1],  # Backward kcat
                'molmass': row['molMass'],
                'genes': [gene_id],
                'complex_with': None
            }
        #else:
            # rxn2protein[rxn_id][enzyme_id]['genes'].append(gene_id)

    #if no enzyme info is found, add dummy enzyme with median kcat and molmass
    for rxn in model.reactions:
        if rxn.id not in rxn2protein.keys() and 'EX'.lower() not in rxn.id.lower() and 'BIOMASS' not in rxn.id and rxn.genes:
            rev = _check_reaction_reversibility(rxn)
            if rev == 0:
                kcat_dict = {'f': 22}
            elif rev ==1:
                kcat_dict = {'b': 22}
            else:
                kcat_dict = {'f': 22,'b': 22}
            # no enzyme information found
            # print('No enzyme information found for reaction: ' + rxn.id)
            enzyme_id = 'Enzyme_'+rxn.id
            rxn2protein[rxn.id] = {enzyme_id:{
                **kcat_dict,
                'molmass': 3.947778784340140e04}}
            #add geneinfo for unknown enzymes
            gpr_info = parse_gpr_information_for_protein2genes(rxn.gpr, enzyme_id)
            protein2gene[enzyme_id] = gpr_info

            for gene in [gene for sublist in gpr_info for gene in sublist]:
                if gene not in gene2transcript.keys():
                    gene2transcript[gene] = {'id': f'mRNA_{gene}', 'length': 750.0}

    return rxn2protein, protein2gene, gene2transcript

def _check_rxn_identifier_format(rxn_id:str) -> str:
    if 'pp' in rxn_id:
        idx_new = rxn_id.replace('pp', '')
    elif 'ex' in rxn_id:
        idx_new = rxn_id.replace('ex', '')
    else:
        idx_new = rxn_id
    return idx_new

def _get_fwd_bckw_kcat(rxn_id: str, kcat:float, model:PAModel) -> Union[list, None]:
    #skip exchange and biomass reaction
    if 'EX' in rxn_id or 'BIOMASS' in rxn_id:
        return None

    # Extract the base identifier without any suffixes
    base_id = rxn_id.split('_')[0]

    # Iterate over each identifier in the input
    if base_id in model.reactions:
        # Determine the form of the identifier
        if rxn_id.endswith('_f'):
            kcat_fwd = kcat
            kcat_rev = 0
        elif rxn_id.endswith('_b'):
            kcat_fwd = 0
            kcat_rev = kcat
        # identifier has no suffix
        elif base_id == rxn_id:
            kcat_fwd = kcat
            kcat_rev = kcat
        else:
            return None
    elif rxn_id in model.reactions:
        kcat_fwd = kcat
        kcat_rev = kcat
    else:
        return None
    return [kcat_fwd, kcat_rev]

def _check_reaction_reversibility(reaction):
    if reaction.lower_bound >= 0:
        # irreversible reaction (forward direction)
        rev = 0
    elif reaction.upper_bound <= 0:
        # irreversible reaction (reverse direction)
        rev = 1
    else:
        rev = 2
        # reversible r
    return rev

def get_kcat_info_for_reaction(enzyme_db: pd.DataFrame, rxn: cobra.Reaction) -> dict:
    rev = 0  # bool to indicate reversibility
    enzyme_db_for_reaction = enzyme_db[enzyme_db['rxnID'] == rxn.id]
    kcat_dict = {}
    if rxn.lower_bound >= 0:
        # irreversible reaction (forward direction)
        rev = 0
        if rxn.id in enzyme_db.rxnID:
            kcat_dict = {'f': enzyme_db_for_reaction.loc['kcat']}
    elif rxn.upper_bound <= 0:
        # irreversible reaction (reverse direction)
        rev = 1
        rxn_id = rxn.id + '_b'
        if rxn_id in enzyme_db.rxnID:
            kcat_dict = {'b': enzyme_db.loc[rxn_id, 'kcat']}
    else:
        rev = 2
        # reversible reaction
        rxn_id_f = rxn.id + '_f'
        rxn_id_b = rxn.id + '_b'
        if rxn_id_f in enzyme_db.rxnID and rxn_id_b in enzyme_db.rxnID:
            kcat_dict = {'f': enzyme_db[enzyme_db['rxnID'] == rxn_id_f]['kcat'],
                             'b': enzyme_db[enzyme_db['rxnID'] == rxn_id_b]['kcat']}

        else:
            # try if only forward reaction is in database
            kcat_dict = {'f': enzyme_db_for_reaction['kcat'],
                             'b': enzyme_db_for_reaction[
                                      'kcat'] / 2}  # deduce backwards kcat from forward value

    return kcat_dict, rev


def parse_gpr_information_for_protein2genes(gpr_info:str, enzyme_id):
    #filter out nan entries
    if isinstance(gpr_info, str):
        return parse_and_or_gpr_relations_from_string(gpr_info)

        # nested_list = ast.literal_eval(gpr_info)
        # gpr_list += split_list(nested_list, seperator=' or ')
    else:
        return [['gene_'+enzyme_id]]

def parse_and_or_gpr_relations_from_string(gpr_info:str):
    gpr_list = []
    for gene_relation in gpr_info.split(' or '):
        sublist = [gene.strip("(')") for gene in gene_relation.split(' and ')]
        gpr_list+=[sublist]
    return gpr_list
#
# def split_list(original_list, seperator):
#     # Create a list to hold the split lists
#     split_lists = []
#     part_list = []
#     # Iterate over the elements of the original list
#     for i, sublist in enumerate(original_list):
#         if isinstance(sublist, list):
#             for string in sublist:
#                 part_list, split_lists = parse_split_list_string_entry(string, seperator, part_list, split_lists)
#         else:
#             part_list, split_lists = parse_split_list_string_entry(sublist, seperator, part_list, split_lists)
#     split_lists += [part_list]
#     return split_lists
#
# def parse_split_list_string_entry(string, seperator, part_list, split_lists):
#     if seperator in string:
#         split_items = string.split(seperator)
#         part_list += [split_items[0].strip("(')")]
#         split_lists += [part_list]
#         part_list = [split_items[1].strip("(')")]
#     else:
#         part_list += [string.strip("(')")]
#     return part_list, split_lists
