from Scripts.mcpam_simulations_analysis import (run_simulations_pam_mcpam_w_different_areas,
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
    pam_info_path = 'Results/PAM_parametrizer/Enzymatic_files/2026_03_11/proteinAllocationModel_EnzymaticData_mciML1515_biomass_weight_8_5.xlsx'
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

    run_simulation_pam_mcpam(models=models)


    










