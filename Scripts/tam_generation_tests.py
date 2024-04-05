import os
import pandas as pd
import numpy as np
import cobra
from typing import Union
import ast

from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.PAModelpy.configuration import Config

from src.TAModelpy import TAModel, ActivemRNASector

AVERAGE_MRNA_LENGTH = 750 #nt

def set_up_ecolicore_tam_simple(total_protein:bool = True, total_mrna = True, active_enzymes: bool = True,
                         translational_enzymes:bool = True, unused_enzymes:bool = True,
                         sensitivity =True):
    # Setting the relative paths
    DATA_DIR = os.path.join('Data')
    MODEL_DIR = os.path.join('Models')
    TAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'TAModel','2024-02-27_gene_enzyme_reaction_relation_Ecoli.xlsx')


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
        rxn2protein, protein2gene, gene2transcript = parse_reaction2protein2gene2transcript_simple(enzyme_db, model)

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

def set_up_ecolicore_tam_specified_transcripts(num_transcripts:int = 0, total_protein:bool = True, total_mrna = True, active_enzymes: bool = True,
                         translational_enzymes:bool = True, unused_enzymes:bool = True,
                         sensitivity =True):
    # Setting the relative paths
    DATA_DIR = os.path.join('Data')
    MODEL_DIR = os.path.join('Models')
    TAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'TAModel','2024-02-27_gene_enzyme_reaction_relation_Ecoli.xlsx')


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
        rxn2protein, protein2gene, gene2transcript = parse_reaction2protein2gene2transcript(enzyme_db, model, num_transcripts)

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


def parse_reaction2protein2gene2transcript_simple(enzyme_db: pd.DataFrame, model:cobra.Model) -> dict:
    # Initialize dictionaries
    rxn2protein = {}
    protein2gene = {}
    gene2transcript = {}

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
            gene_id = 'gene_'+enzyme_id
            mrna_length = AVERAGE_MRNA_LENGTH


            # Create gene-protein-reaction associations
            if enzyme_id not in protein2gene:
                protein2gene[enzyme_id] = [[gene_id]]

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
        else:
            rxn2protein[rxn_id][enzyme_id]['genes'].append(gene_id)

    #if no enzyme info is found, add dummy enzyme with median kcat and molmass
    for rxn in model.reactions:
        if rxn.id not in rxn2protein.keys() and 'EX' not in rxn.id and 'BIOMASS' not in rxn.id and rxn.genes:
            rev = _check_reaction_reversibility(rxn)
            if rev == 0:
                kcat_dict = {'f': 22}
            elif rev ==1:
                kcat_dict = {'b': 22}
            else:
                kcat_dict = {'f': 22,'b': 22}
            # no enzyme information found
            print('No enzyme information found for reaction: ' + rxn.id)
            enzyme_id = 'Enzyme_'+rxn.id
            gene_id = 'gene_'+enzyme_id
            rxn2protein[rxn.id] = {enzyme_id:{
                **kcat_dict,
                'molmass': 3.947778784340140e04}}
            #add geneinfo for unknown enzymes
            protein2gene[enzyme_id] =[[gene_id]]

            gene2transcript[gene_id] = {'id': f'mRNA_{gene_id}', 'length': 750.0}

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

def add_transcript_or_relation(tam):
    index = np.random.randint(low=0, high=len(tam.enzymes))
    enzyme_ut = tam.enzymes[index]
    genes_ut = ['test1', 'test2']
    gene_length = [10, 100]
    # Act
    enzyme_ut.add_genes(genes_ut, gene_length)
    return tam

def add_transcript_and_relation(tam):
    index = np.random.randint(low=0,high=len(tam.enzymes))
    enzyme_ut = tam.enzymes[index]
    genes_ut = ['test3', 'test4']
    gene_length = [10, 100]

    # Act
    enzyme_ut.add_genes(genes_ut, gene_length, relation='AND')
    print(enzyme_ut)
    return tam


def parse_reaction2protein2gene2transcript(enzyme_db: pd.DataFrame, model:cobra.Model, num_transcripts: int = 1000) -> dict:
    # Initialize dictionaries
    rxn2protein = {}
    protein2gene = {}
    gene2transcript = {}

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
            if not isinstance(gene_id, str): gene_id = 'gene_'+enzyme_id
            mrna_length = row['mrna_length']

            #get the gene-protein-reaction-associations
            gpr = parse_gpr_information_for_protein2genes(row['gpr'], enzyme_id)

            # Create gene-protein-reaction associations
            if enzyme_id not in protein2gene:
                protein2gene[enzyme_id] = gpr

            # Create gene2transcript dictionary
            if gene_id not in gene2transcript.keys() and len(gene2transcript) < num_transcripts:
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
        else:
            rxn2protein[rxn_id][enzyme_id]['genes'].append(gene_id)

    #if no enzyme info is found, add dummy enzyme with median kcat and molmass
    for rxn in model.reactions:
        if rxn.id not in rxn2protein.keys() and 'EX' not in rxn.id and 'BIOMASS' not in rxn.id and rxn.genes:
            rev = _check_reaction_reversibility(rxn)
            if rev == 0:
                kcat_dict = {'f': 22}
            elif rev ==1:
                kcat_dict = {'b': 22}
            else:
                kcat_dict = {'f': 22,'b': 22}
            # no enzyme information found
            print('No enzyme information found for reaction: ' + rxn.id)
            enzyme_id = 'Enzyme_'+rxn.id
            rxn2protein[rxn.id] = {enzyme_id:{
                **kcat_dict,
                'molmass': 3.947778784340140e04}}
            #add geneinfo for unknown enzymes
            gpr_info = parse_gpr_information_for_protein2genes(rxn.gpr, enzyme_id)
            protein2gene[enzyme_id] = gpr_info

            for gene in [gene for sublist in gpr_info for gene in sublist]:
                if gene not in gene2transcript.keys() and len(gene2transcript) < num_transcripts:
                    gene2transcript[gene] = {'id': f'mRNA_{gene}', 'length': 750.0}

    return rxn2protein, protein2gene, gene2transcript

def parse_gpr_information_for_protein2genes(gpr_info:str, enzyme_id):
    #filter out nan entries
    if isinstance(gpr_info, str):
        nested_list = ast.literal_eval(gpr_info)
        # Extracting the inner lists and removing parentheses
        gpr_list = [[[item.strip("(')")] + sublist[1:] for item in sublist[0].split(" or ")] for sublist in nested_list]
        return gpr_list[0]
    else:
        return [['gene_'+enzyme_id]]


if __name__ == '__main__':
    tam = set_up_ecolicore_tam_simple(total_mrna=True)
    tam.optimize()
    print('Predicted growth rate: ',tam.objective.value)
    # OR
    print('\nAdding a single OR relationship: test1 OR test2 to a random enzyme in the model')
    tam = add_transcript_or_relation(tam)
    tam.optimize()
    print('Predicted growth rate: ',tam.objective.value)

    #AND
    print('\nAdding a single AND relationship: test3 AND test4 to a random enzyme in the model')
    tam = add_transcript_and_relation(tam)
    tam.optimize()
    # for constr in tam.transcripts.get_by_id('mRNA_test3_test4')._constraints.values():
    #     print(constr, 'primal: ', constr.primal, 'dual: ', constr.dual)
    print('Predicted growth rate: ',tam.objective.value)

    print('\n Building the wildtype ecolicore tam with a single transcript added to it')
    tam = set_up_ecolicore_tam_specified_transcripts(num_transcripts = 100, total_mrna=True)
    tam.optimize()
    print('Predicted growth rate: ',tam.objective.value)
    # trans = tam.transcripts.get_by_id('mRNA_b0720')
    # for key, value in tam.constraints.items():
    #     if 'mRNA' in key: print(value)
    # for constr in tam.transcripts.get_by_id('mRNA_b0720')._constraints.values():
    #     print(constr, 'primal: ', constr.primal, 'dual: ', constr.dual)
    # print(tam.transcripts)