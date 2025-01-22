from Scripts.mcpam_simulations_analysis import (run_pam_mcpam_core_with_optimized_kcats,
                                                 run_simulation_pam_mcpam,
                                                 run_simulations_pam_mcpam_w_different_areas,
                                                 set_up_ecolicore_pam,
                                                 set_up_ecolicore_mcpam,
                                                 set_up_ecolicore_mcpam_new_surface_parameter,
                                                 compare_mu_for_different_sensitivities_ecolicore_pam,
                                                 perform_single_gene_ko_for_all_genes,
                                                 perform_and_plot_single_KO)
from Scripts.mcpam_generation_uniprot_id import set_up_ecoli_mcpam, set_up_ecoli_pam
from Scripts.mcpam_toy_generation import build_toy_model
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
import matplotlib.pyplot as plt; plt.rcdefaults()

if __name__ == "__main__":
    # Build full scale pam new parsing
    # pam_info_path = 'Data/proteinAllocationModel_mciML1515_EnzymaticData_241209_multi.xlsx'
    # pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=False)
    # mcpam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)
    #
    # models = [pam, mcpam]
    # run_simulation_pam_mcpam(models, type="full scale")


    # Build core pam
    # pam = set_up_ecolicore_pam(sensitivity=False)
    # mcpam = set_up_ecolicore_mcpam_new_surface_parameter(sensitivity=False)
    # models = [pam, mcpam]
    # run_simulation_pam_mcpam(models, type="core")

    # Build toy pam
    toy_pam = build_toy_model(membrane_sector=True)
    toy_pam.objective = "R11"
    toy_pam.optimize()
    print(toy_pam.objective.value)











