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

    ## Build full scale pam and change the enzyme sectors accordingly (based on script from Tobias. A)

    # Load iML1515 PAM and mcPAM
    pam_info_path = 'Data/proteinAllocationModel_EnzymaticData_iML1515_10.xlsx'
    model_path = 'Models/iML1515.xml'
    pam = set_up_pam(pam_info_file=pam_info_path,
                    model=model_path,
                    sensitivity=False,
                    membrane_sector=False,                  
                    )
    mcpam = set_up_pam(pam_info_file=pam_info_path, 
                       model=model_path,
                       sensitivity=False, 
                       membrane_sector=True,
                       separate_memprot_from_tpc=True,
                       total_protein=0.1935
                       )
    models = [pam, mcpam]

    for tot_prot in [0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20]:
        pam.change_total_protein_constraint(tot_prot)
        mcpam.change_total_protein_constraint(tot_prot)
        run_simulations_pam_mcpam_w_different_areas(models, type='full scale', max_area_list=[0.03, 0.04, 0.15, 0.50, 1])


    










