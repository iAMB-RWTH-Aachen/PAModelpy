from Scripts.mcpam_simulations_analysis import (run_pam_mcpam_core_with_optimized_kcats,
                                                 run_simulation_pam_mcpam,
                                                 run_simulations_pam_mcpam_w_different_areas,
                                                 set_up_ecolicore_pam,
                                                 set_up_ecolicore_mcpam,
                                                 set_up_ecolicore_mcpam_new_surface_parameter,
                                                 compare_mu_for_different_sensitivities_ecolicore_pam,
                                                 perform_single_gene_ko_for_all_genes,
                                                 perform_and_plot_single_KO,
                                                 get_memprot_data_in_mcpam)
from Scripts.create_pamodel_from_diagnostics_file import (create_pamodel_from_diagnostics_file,
                                                          change_translational_sector_with_config_dict,
                                                          _set_up_pamodel_for_simulations,
                                                          change_memprot_kcats
                                                          )
from Scripts.mcpam_generation_uniprot_id import set_up_ecoli_mcpam, set_up_ecoli_pam
from Scripts.mcpam_toy_generation import build_toy_model
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
import matplotlib.pyplot as plt; plt.rcdefaults()
import pandas as pd
import re
import os

if __name__ == "__main__":
    # Build full scale old parsing, TS enzymatic dataset
    # pam = set_up_ecoli_pam(sensitivity=False)
    # mcpam = set_up_ecoli_mcpam(sensitivity=False)
    # models = [pam, mcpam]
    #
    # run_simulations_pam_mcpam_w_different_areas(models, type='full scale')


    #Build full scale pam new parsing
    # pam_info_path = 'Results/PAM_parametrizer/Files/2025_02_20/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
    # sheet_name = "diagnostics_1"
    # pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=False)
    # mcpam = set_up_pam(pam_info_file=pam_info_path, sheet_name= sheet_name, sensitivity=False, membrane_sector=True)
    # models = [pam, mcpam]
    # membrane_proteins = mcpam.sectors.get_by_id('MembraneSector').membrane_proteins
    # for enzyme_id, (list, alpha_number) in membrane_proteins.items():
    #     if alpha_number >= 28:
    #         print(enzyme_id, list, alpha_number)
    # run_simulations_pam_mcpam_w_different_areas(models, type="full scale")

    # Build full scale pam from diagnostics file
    ## Define necessary paths/sheet names
    diagnostics_data_path = 'Results/PAM_parametrizer/Files/2025_03_11/pam_parametrizer_diagnostics_mciML1515_1.xlsx'
    pam_info_path = 'Results/PAM_parametrizer/Enzymatic_files/2025_05_14/proteinAllocationModel_EnzymaticData_iML1515_1.xlsx'
    sheet_name = 'Best_Individuals'

    pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=False)
    # _set_up_pamodel_for_simulations(pam, 'EX_glc__D_e', transl_sector_config=True)
    # pam = create_pamodel_from_diagnostics_file(diagnostics_data_path, pam, sheet_name)


    mcpam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)
    # _set_up_pamodel_for_simulations(mcpam, 'EX_glc__D_e', transl_sector_config=True)
    # mcpam = create_pamodel_from_diagnostics_file(diagnostics_data_path, mcpam, sheet_name)
    models = [pam, mcpam]

    # Changing the kcats using excel sheet
    # for model in models:
    #     memprot_file_path = 'Results/PAM_parametrizer/Files/2025_03_11/memprot_data.xlsx'
    #     memprot_sheet_name = 'diagnostics_2'
    #     model = change_memprot_kcats(memprot_file_path, model, memprot_sheet_name)

    mcpam.objective = mcpam.configuration.BIOMASS_REACTION
    solution = mcpam.optimize()
    mcpam.sectors.get_by_id('MembraneSector').calculate_occupied_membrane(mcpam)
    print(solution.fluxes.ACt2rpp)
    print(solution.fluxes.ACt4pp)
    print(solution.fluxes.ACtex)

    # # Get a df for all proteins and their occupancy in the membrane
    # file_name = os.path.basename(pam_info_path)  # 'proteinAllocationModel_EnzymaticData_iML1515_10.xlsx'
    # file_stem = os.path.splitext(file_name)[0]  # 'proteinAllocationModel_EnzymaticData_iML1515_10'
    # match = re.search(r'_(\d+)$', file_stem)
    # number = match.group(1) if match else None
    # prot_occupancy_df = mcpam.sectors.get_by_id('MembraneSector').calculate_occupied_membrane(mcpam, get_df=True)

    # # Write excel datasheet
    # data_path = os.path.join('Results/PAM_parametrizer/Enzymatic_files/2025_05_14/protein_occupancy_data.xlsx')
    # with pd.ExcelWriter(data_path, engine='openpyxl', mode='a') as writer:
    #     # Write the new DataFrame to a new sheet
    #     prot_occupancy_df.to_excel(writer, sheet_name=f'enzymatic_file_{number}', index=True)

    # # run_simulation_pam_mcpam(models, type='full scale')
    # run_simulations_pam_mcpam_w_different_areas(models, type="full scale")

    # # Build core pam
    # pam = set_up_ecolicore_pam(sensitivity=False)
    # mcpam = set_up_ecolicore_mcpam_new_surface_parameter(sensitivity=False)
    # # # Get a df for all proteins and their occupancy in the membrane
    # mcpam.objective = mcpam.configuration.BIOMASS_REACTION
    # mcpam.optimize()
    # mcpam.sectors.get_by_id('MembraneSector').calculate_occupied_membrane(mcpam, get_df=True)
    # models = [pam, mcpam]
    # run_simulations_pam_mcpam_w_different_areas(models, type="core")

    # Build core pam new parse
    # pam = set_up_core_pam(sensitivity=False, membrane_sector=False, pam_info_file='Data/proteinAllocationModel_mc-core_EnzymaticData_241209_multi.xlsx')
    # mcpam = set_up_core_pam(sensitivity=False, membrane_sector=True, pam_info_file='Data/proteinAllocationModel_mc-core_EnzymaticData_241209_multi.xlsx')
    # models = [pam, mcpam]
    # run_simulation_pam_mcpam(models, type="core")

    # Build toy pam
    # toy_pam = build_toy_model(membrane_sector=True)
    # toy_pam.objective = "R11"
    # toy_pam.optimize()
    # print(toy_pam.objective.value)











