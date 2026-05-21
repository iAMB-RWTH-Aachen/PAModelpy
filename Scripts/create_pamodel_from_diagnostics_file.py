import pandas as pd
from cobra import Model
from typing import Tuple, Literal
import re
import math
import numpy as np
from typing import Union
from src.PAModelpy.PAModel import PAModel, ActiveEnzymeSector, MembraneSector
from src.PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein, _order_enzyme_complex_id
from Scripts.mcpam_simulations_analysis import (
    run_simulation_pam_mcpam, 
    run_simulations_pam_mcpam_w_different_areas,
    run_simulations_pam_mcpam_w_different_tpcs
)

DEFAULT_MOLMASS = 39959.4825 #kDa
DEFAULT_KCAT = 11 #s-1

def _set_up_pamodel_for_simulations(pamodel:PAModel,
                                   substrate_id: str,
                                   transl_sector_config:Union[bool, dict[str, float]]) -> None:
    if not isinstance(transl_sector_config, dict) and transl_sector_config:
        transl_sector_config = {'slope': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                                'intercept': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_0[0]}

    if transl_sector_config is not False:
        change_translational_sector_with_config_dict(pamodel=pamodel,
                                                     transl_sector_config = transl_sector_config,
                                                     substrate_uptake_id = substrate_id)

def change_translational_sector_with_config_dict(pamodel:PAModel,
                                                 transl_sector_config:dict,
                                                 substrate_uptake_id:str) -> None:
    pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].lb = 0 #need to set the lb to 0 to prevent errors in the setter methods

    pamodel.change_sector_parameters(pamodel.sectors.get_by_id('TranslationalProteinSector'),
                                              slope=transl_sector_config['slope'],
                                              intercept=transl_sector_config['intercept'],
                                    lin_rxn_id=substrate_uptake_id
                                     )

def change_prot_kcats(prot_df:pd.DataFrame, model:Union[Model, PAModel])-> Union[Model, PAModel]:
    
    for _, row in prot_df.iterrows():
        if row['Reaction'].startswith('CE_'):
            rxn_id = _extract_reaction_id_from_catalytic_reaction_id(row['Reaction'])
        else:
            rxn_id = row['Reaction']
        enzyme_id = _order_enzyme_complex_id(row['enzyme_id'])
        kcat_dict = {rxn_id: {'f': row['Forward Flux'], 'b': row['Backward Flux']}}
    
        model.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)

    return model

def create_pamodel_from_diagnostics_file(file_path:str,
                                         model: PAModel,
                                         sheet_name: str = 'Best_Individuals',
                                         enzyme_sector_update: bool = True,
                                         other_enzyme_id_pattern: str = r'E[0-9][0-9]*|Enzyme_\S+',
                                         substrate_uptake_id: str = 'EX_glc__D_e'
                                         )-> PAModel:
    """
    Modifies a Protein Allocation Model using information about turnover numbers from a diagnostics file
    (result from PAMparametrizer). If available, also adjusts the sector parameters associated to the substrate.

    Args:
        file_path (str): path to the diagnostics xlsx file. The file should at least have the following columns:
            - run_id: the iteration of the PAMparametrizer
            - rxn_id: the id of the reaction to modify. This can be the catalytic reaction id (CE_<rxn_id>_<enzyme_id>)
            - enzyme_id: the id of the enzyme to modify
            - direction: 'f' or 'b', determines the directionality of the reaction
            - kcat[s-1]: the new kcat value in 1/s
        model (PAModel): the PAM to adjust
        sheet_name (str): name of the sheet with the information about the modifications
        enzyme_sector_update (bool): if the enzyme sectors should be updated according to the parametrization results.
            Defaults to True
        other_enzyme_id_pattern (regex str): regex pattern which matches default enzyme id.
            Normally is E1, E2, etc. or Enzyme_<rxn_id>, but can be specified by the user.
        substrate_uptake_id (str): name of the uptake reaction for the substrate for which the model is built

    Returns:
        PAModel: with adjusted parameters

    """
    best_individual_df = pd.read_excel(file_path, sheet_name=sheet_name)
    for _, group in best_individual_df.groupby('run_id'):
        for _, row in group.iterrows():
            rxn_id = _extract_reaction_id_from_catalytic_reaction_id(row['rxn_id'],
                                                                     default_enzyme_id_pattern = other_enzyme_id_pattern
                                                                     )
            enzyme_id = _order_enzyme_complex_id(row['enzyme_id'],
                                                 other_enzyme_id_pattern = other_enzyme_id_pattern)
            kcat_dict = {rxn_id: {row['direction']: row['kcat[s-1]']}}
            model.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)
            
    if not enzyme_sector_update: return model
    try:
        sector_parameters_df = pd.read_excel(file_path, sheet_name="sector_parameters")
    except:
        return model

    for sector, sector_params in sector_parameters_df.groupby('sector_id'):
        sector_params = sector_params.loc[
            (sector_parameters_df.substrate_uptake_id == substrate_uptake_id)
        ].rename({'substrate_uptake_id': 'lin_rxn_id'},
                 axis=1)[['slope', 'intercept', 'lin_rxn_id']].to_dict('records')[0]
        model.change_sector_parameters(
            sector = model.sectors.get_by_id(sector),
            **sector_params,
            print_change=True
        )
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

def _extract_reaction_id_from_catalytic_reaction_id(input_str: str,
                                                    default_enzyme_id_pattern: str = r'E[0-9][0-9]*|Enzyme_*') -> str:
    # Define the regex pattern for protein IDs, obtained from UniProtKB, 2024-08-07
    # https://www.uniprot.org/help/accession_numbers
    protein_id_pattern = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'

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

if __name__ == '__main__':
    pam_info_file = 'Results/PAM_parametrizer/Diagnostics_files/2026_05_21/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
    model_path = 'Models/iML1515.xml'
    pam = set_up_pam(pam_info_file=pam_info_file,
                    model=model_path,
                    sensitivity=False,
                    membrane_sector=False,                  
                    )
    mcpam = set_up_pam(pam_info_file=pam_info_file, 
                       model=model_path,
                       sensitivity=False, 
                       membrane_sector=True,
                       separate_memprot_from_tpc=True,
                       total_protein=0.241,
                       max_membrane_area=0.5154
                       )

    mcpam = create_pamodel_from_diagnostics_file(file_path='Results/PAM_parametrizer/Diagnostics_files/2026_05_21/pam_parametrizer_diagnostics_mciML1515_1.xlsx',
                                                 model=mcpam,
                                                 sheet_name='Best_Individuals')
    # mcpam.change_sector_parameters(mcpam.unused_enzymes, slope=0.0075, intercept=None, lin_rxn_id='EX_glc__D_e')
    # sector_parameters_df = pd.read_excel('Results/PAM_parametrizer/Diagnostics_files/2026_05_06/pam_parametrizer_diagnostics_mciML1515_2.xlsx', sheet_name="sector_parameters")

    # for sector, sector_params in sector_parameters_df.groupby('sector_id'):
    #     sector_params = sector_params.loc[
    #         (sector_parameters_df.substrate_uptake_id == 'EX_glc__D_e')
    #     ].rename({'substrate_uptake_id': 'lin_rxn_id'},
    #              axis=1)[['slope', 'intercept', 'lin_rxn_id']].to_dict('records')[0]
    #     model.change_sector_parameters(
    #         sector = model.sectors.get_by_id(sector),
    #         **sector_params,
    #         print_change=True
    #     )
    
    models = [pam, mcpam]

    run_simulation_pam_mcpam(models)    
    # run_simulations_pam_mcpam_w_different_areas(models=models, max_area_list=[0.2, 0.3, 0.4, 0.5, 0.6])
    # run_simulations_pam_mcpam_w_different_tpcs(models=models, tpc_list=[0.2, 0.21, 0.22, 0.23, 0.24, 0.258])
    
    mcpam.reactions.get_by_id('EX_glc__D_e').lower_bound = -10
    mcpam.reactions.get_by_id('EX_glc__D_e').upper_bound = -10
    mcpam.optimize()
    print(mcpam.objective.value)
    
    occupied_area, available_area = mcpam.sectors.get_by_id('MembraneSector').calculate_occupied_membrane(mcpam)
    print(available_area, occupied_area)
    print(mcpam.solver.shadow_prices['membrane'])
    print(mcpam.constraints['membrane'].dual)



    