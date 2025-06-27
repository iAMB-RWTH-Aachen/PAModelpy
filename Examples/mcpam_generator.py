from Scripts.mcpam_simulations_analysis import (run_simulations_pam_mcpam_w_different_areas,
                                                run_simulation_pam_mcpam,
                                                get_info_for_proteins,
                                                get_missing_backward_kcats,
                                                fill_missing_backward_kcats)
from Scripts.create_pamodel_from_diagnostics_file import (create_pamodel_from_diagnostics_file,
                                                          change_translational_sector_with_config_dict,
                                                          _set_up_pamodel_for_simulations,
                                                          change_prot_kcats
                                                          )
from Scripts.mcpam_generation_uniprot_id import set_up_ecoli_mcpam, set_up_ecoli_pam
from Scripts.mcpam_toy_generation import build_toy_model
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
import matplotlib.pyplot as plt; plt.rcdefaults()
import pandas as pd
import re
import os

if __name__ == "__main__":

    ## Build full scale pam and change the enzyme sectors accordingly (based on script from Tobias. A)
    

    ### Compare mcPAM using full scale GEM for kcat dataset with missing backward kcats vs. dataset with filled backward kcats.
    # diff_mu = []

    # for i in range(0,10):
    #     number = i + 1
    #     pam_info_path = f'Results/PAM_parametrizer/Enzymatic_files/2025_05_14/proteinAllocationModel_EnzymaticData_iML1515_{number}.xlsx'
    #     mcpam_missing_kcat_b = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)
    #     mcpam_filled_kcat_b = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)

    #     missing_backward_kcats = get_missing_backward_kcats(mcpam_missing_kcat_b)
    #     filled_backward_kcats = fill_missing_backward_kcats(missing_backward_kcats)

    #     mcpam_filled_kcat_b = change_prot_kcats(prot_df=filled_backward_kcats, model=mcpam_filled_kcat_b)
        
    #     models = [mcpam_missing_kcat_b, mcpam_filled_kcat_b]

    #     fig = run_simulation_pam_mcpam(models=models)

    #     plt.savefig(f"Results/PAM_parametrizer/Analysis/Flux_simulation_mcpams_comparison_missing_and_filled_kcats_{number}.png", dpi=300)


    ### Build full scale pam from diagnostics file
    ## Define necessary paths/sheet names
    # diagnostics_data_path = 'Results/PAM_parametrizer/Files/2025_03_11/pam_parametrizer_diagnostics_mciML1515_1.xlsx'
    # pam_info_path = 'Results/PAM_parametrizer/Enzymatic_files/2025_05_14/proteinAllocationModel_EnzymaticData_iML1515_10.xlsx'
    # sheet_name = 'Best_Individuals'

    # pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=False)
    # _set_up_pamodel_for_simulations(pam, 'EX_glc__D_e', transl_sector_config=True) # changing the translational sector
    # pam = create_pamodel_from_diagnostics_file(diagnostics_data_path, pam, sheet_name)


    # mcpam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)
    # mcpam.optimize()
    # _set_up_pamodel_for_simulations(mcpam, 'EX_glc__D_e', transl_sector_config=True) # changing the translational sector
    # mcpam = create_pamodel_from_diagnostics_file(diagnostics_data_path, mcpam, sheet_name)
    # models = [pam, mcpam]

    # get_info_for_proteins(mcpam=model,
    #                       pam_info_path=pam_info_path,
    #                       protein_info_path="Results/PAM_parametrizer/Enzymatic_files/2025_05_14/protein_occupancy_data.xlsx")

    # change_set_of_kcats_using_excel_sheet(models=models, 
    #                                       prot_file_path="Results/PAM_parametrizer/Enzymatic_files/2025_05_14/protein_occupancy_data.xlsx",
    #                                       sheet="enzymatic_file_10")

    # # run_simulation_pam_mcpam(models, type='full scale')
    # run_simulations_pam_mcpam_w_different_areas(models, type="full scale")










