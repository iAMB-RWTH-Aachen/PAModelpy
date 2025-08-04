from Scripts.mcpam_simulations_analysis import (change_set_of_kcats_using_excel_sheet,
                                                run_simulations_pam_mcpam_w_different_areas,
                                                run_simulation_pam_mcpam,
                                                get_info_for_proteins,
                                                get_missing_backward_kcats,
                                                fill_missing_backward_kcats)
from Scripts.create_pamodel_from_diagnostics_file import (create_pamodel_from_diagnostics_file,
                                                          change_translational_sector_with_config_dict,
                                                          _set_up_pamodel_for_simulations
                                                          )
from Scripts.mcpam_generation_uniprot_id import set_up_ecoli_mcpam, set_up_ecoli_pam
from Scripts.mcpam_toy_generation import build_toy_model
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
import matplotlib.pyplot as plt; plt.rcdefaults()
import pandas as pd
import re
import os

if __name__ == "__main__":

    ## Build full scale pam from diagnostics file
    # Define necessary paths/sheet names
    diagnostics_data_path = 'Results/PAM_parametrizer/Diagnostics_files/2025_07_24/pam_parametrizer_diagnostics_mciML1515_2.xlsx'
    pam_info_path = 'Results/PAM_parametrizer/Enzymatic_files/2025_07_24/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
    sheet_name = 'Best_Individuals'

    pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=False)
    _set_up_pamodel_for_simulations(pam, 'EX_glc__D_e', transl_sector_config=True) # changing the translational sector
    pam = create_pamodel_from_diagnostics_file(diagnostics_data_path, pam, sheet_name)


    mcpam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)
    _set_up_pamodel_for_simulations(mcpam, 'EX_glc__D_e', transl_sector_config=True) # changing the translational sector
    mcpam = create_pamodel_from_diagnostics_file(diagnostics_data_path, mcpam, sheet_name)
    models = [pam, mcpam]

    # for model in models:
    #     ue_sector = model.sectors.get_by_id('UnusedEnzymeSector')
    #     te_sector = model.sectors.get_by_id('TranslationalProteinSector')
    #     # Change unused enzyme sector parameters
    #     model.change_sector_parameters(sector = ue_sector,
    #                                 slope = 0.013452, #in this case: g_p*h/(g_cdw*mmol_glc) 0.01307
    #                                 intercept=0.172231, # g_p/g_cdw
    #                                 lin_rxn_id= 'EX_glc__D_e', # the reaction that is used to calculate the slope
    #                                 print_change = True #do you want to see the change? False by default
    #                                 )
    #     # Change translational enzyme sector parameters
    #     model.change_sector_parameters(sector = te_sector,
    #                                 slope = -0.00434, #in this case: g_p*h/(g_cdw*mmol_glc)
    #                                 intercept=0.046137, # g_p/g_cdw
    #                                 lin_rxn_id= 'EX_glc__D_e', # the reaction that is used to calculate the slope, EX_glc__D_e
    #                                 print_change = True #do you want to see the change? False by default
    #                                 )

    # get_info_for_proteins(mcpam=model,
    #                       pam_info_path=pam_info_path,
    #                       protein_info_path="Results/PAM_parametrizer/Enzymatic_files/2025_05_14/protein_occupancy_data.xlsx")

    # change_set_of_kcats_using_excel_sheet(model=model, 
    #                     prot_file_path="Results/From_kcat_dataset_20250627/kcats_sensitive_membrane_enzymes_manually_curated.xlsx",
    #                     sheet="bulky_ME")

    # run_simulation_pam_mcpam(models, type='full scale')
    run_simulations_pam_mcpam_w_different_areas(models, type="full scale")











