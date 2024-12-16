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

    pam = set_up_ecoli_pam(sensitivity=False)
    mcpam = set_up_ecoli_mcpam(sensitivity=False)

    # pam, mcpam = run_pam_mcpam_core_with_optimized_kcats(sensitivity=True, type='full scale', enzyme_sets_name="enzyme_sets.json")
    models = [pam, mcpam]
    run_simulations_pam_mcpam(models, type="full scale")

    for protein, kcats in mcpam.sectors.get_by_id("MembraneSector").membrane_proteins.items():
        print(protein, kcats)










