from Scripts.mcpam_simulations_analysis import (run_pam_mcpam_core_with_optimized_kcats,
                                                 run_simulations_pam_mcpam,
                                                 set_up_ecolicore_pam,
                                                 set_up_ecolicore_mcpam,
                                                 compare_mu_for_different_sensitivities_ecolicore_pam,
                                                 perform_single_gene_ko_for_all_genes,
                                                 perform_and_plot_single_KO)
from Scripts.mcpam_generation_uniprot_id import set_up_ecoli_mcpam, set_up_ecoli_pam
import matplotlib.pyplot as plt; plt.rcdefaults()

if __name__ == "__main__":

    mcpam = set_up_ecolicore_mcpam(sensitivity=False)
    pam = set_up_ecolicore_pam(sensitivity=False)
    pam.optimize()
    mcpam.optimize()
    print("PAM objective value: ", pam.objective.value)
    print("mcPAM objective value: ", mcpam.objective.value)

    models = [pam, mcpam]
    run_simulations_pam_mcpam(models, type="core")








