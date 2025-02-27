import pandas as pd
from typing import Tuple, Literal
import re

from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein, _order_enzyme_complex_id

DEFAULT_MOLMASS = 39959.4825 #kDa
DEFAULT_KCAT = 11 #s-1

def create_pamodel_from_diagnostics_file(file_path:str, model: PAModel, sheet_name: str)-> PAModel:
    best_individual_df = pd.read_excel(file_path, sheet_name=sheet_name)
    for _, group in best_individual_df.groupby('run_id'):
        for _, row in group.iterrows():
            rxn_id = _extract_reaction_id_from_catalytic_reaction_id(row['rxn_id'])
            enzyme_id = _order_enzyme_complex_id(row['enzyme_id'])
            kcat_dict = {rxn_id: {row['direction']: row['kcat[s-1]']}}
            model.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)
    return model

def get_rxn2kcat_protein2gene_dict(param_file_path:str, model_file_path: str
                                   ) -> Tuple[
    dict[str, dict[str,dict[Literal['f', 'b', 'molmass', 'protein_reaction_association'], float]]],
    dict[str,str]]:
    # create enzyme objects for each gene-associated reaction
    pam = set_up_pam(param_file_path, model_file_path, sensitivity=False)
    enzyme_db = pd.read_excel(param_file_path, sheet_name='ActiveEnzymes').iloc[:, 1:]
    rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, pam)

    ae_sector = pam.sectors.ActiveEnzymeSector
    new_rxn2prot = ae_sector.rxn2protein.copy()
    for rxn, enz_dict in ae_sector.rxn2protein.items():
        if rxn[:2]=='CE': continue

        for enzyme_id, enzyme_dict in enz_dict.items():
            protein_reaction = enzyme_dict['protein_reaction_association']
            if not ae_sector._enzyme_is_enzyme_complex(protein_reaction, enzyme_id): continue

            for pr in protein_reaction:
                if not len(pr) > 1: continue

                enzyme_complex_id = '_'.join(pr)
                new_rxn2prot[rxn] = {**new_rxn2prot[rxn],
                                     **{enzyme_complex_id: enzyme_dict}}
    return new_rxn2prot, protein2gene

def _extract_reaction_id_from_catalytic_reaction_id(input_str: str) -> str:
    # Define the regex pattern for protein IDs, obtained from UniProtKB, 2024-08-07
    # https://www.uniprot.org/help/accession_numbers
    protein_id_pattern = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'
    default_enzyme_id_pattern:str = r'E[0-9][0-9]*'

    # Remove the 'CE_' prefix if it exists
    if input_str.startswith('CE_'):
        input_str = input_str[3:]

    # Define the regex pattern to match protein IDs
    protein_id_regex = re.compile(r'_' + protein_id_pattern + r'|' + r'_' + default_enzyme_id_pattern)
     # split off all protein ids from the reaction
    reaction_id = protein_id_regex.split(input_str)[0]
    # Remove any trailing or leading underscores that might remain
    return reaction_id.strip('_')

def _get_rxn2kcat_as_series(rxn2kcat: dict[str, dict],
                                   name: str):

    kcats = {}
    for rxn,enz_dict in rxn2kcat.items():
        for enz, kcat_dict in enz_dict.items():
            for direction, kcat in kcat_dict.items():
                if len(direction) == 1:
                    kcats[f"{rxn}_{enz}_{direction}"] = kcat
    return pd.Series(kcats, name = name)



