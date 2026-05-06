from Scripts.mcpam_simulations_analysis import *
from Scripts.create_pamodel_from_diagnostics_file import (create_pamodel_from_diagnostics_file,
                                                          change_translational_sector_with_config_dict,
                                                          _set_up_pamodel_for_simulations
                                                          )
from Scripts.mcpam_toy_generation import build_toy_model
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
from Scripts.mcpam_generation_uniprot_id import set_up_ecolicore_mcpam_new_surface_parameter, set_up_ecolicore_pam
import matplotlib.pyplot as plt; plt.rcdefaults()
import pandas as pd
import re
import os
from cobra.io import load_json_model

if __name__ == "__main__":

    ## Build full scale pam and change the enzyme sectors accordingly (based on script from Tobias. A)

    # Load iML1515 PAM and mcPAM
    pam_info_path = 'Data/proteinAllocationModel_mc-core_EnzymaticData_241209_multi.xlsx'
    model_path = 'Models/e_coli_core.json'

    core_gem = load_json_model(model_path)
    pam = set_up_core_pam(
                    pam_info_file=pam_info_path,
                    model=model_path,
                    sensitivity=False,
                    membrane_sector=False
                    )
    mcpam = set_up_core_pam(
                        pam_info_file=pam_info_path, 
                       model=model_path,
                       sensitivity=False,
                       membrane_sector=True
                       )
    pam_mcpam = [pam, mcpam]
    models = [core_gem, pam, mcpam]

    # Run simulation for both PAM and mcPAM with the changed sector parameters
    run_simulations_pam_mcpam_w_different_areas(pam_mcpam, type='core')











