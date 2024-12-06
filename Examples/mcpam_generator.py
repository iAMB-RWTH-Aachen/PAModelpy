from Scripts.mcpam_simulations_analysis import (run_pam_mcpam_core_with_optimized_kcats,
                                                 run_simulations_pam_mcpam_core,
                                                 set_up_ecolicore_pam,
                                                 set_up_ecolicore_mcpam,
                                                 compare_mu_for_different_sensitivities_ecolicore_pam,
                                                 perform_single_gene_ko_for_all_genes,
                                                 perform_and_plot_single_KO)
from cobra import Model, Reaction, Metabolite
import pandas as pd
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    mcpam_core = set_up_ecolicore_mcpam(sensitivity=True)
    mcpam_core.optimize()
    print(mcpam_core.objective.value)
    # pam_core = set_up_ecolicore_pam(sensitivity=True)
    #
    # models = [pam_core, mcpam_core]
    #
    # run_simulations_pam_mcpam_core(models)









