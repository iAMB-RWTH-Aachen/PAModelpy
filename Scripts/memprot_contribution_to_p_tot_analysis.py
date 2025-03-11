from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import seaborn as sns

import pandas as pd
import os
import re
import sys
import numpy as np

sys.path.append('C:\\Users\\claud\\Documents\\iamb-student-folders\\iamb-folder-template\\mcPAM_package')

from src.PAModelpy.configuration import Config

if os.path.split(os.getcwd())[1] == 'Figures':
    os.chdir(os.path.split(os.getcwd())[0])
from src.PAModelpy.utils.pam_generation import set_up_pam
from Scripts.create_pamodel_from_diagnostics_file import (create_pamodel_from_diagnostics_file,
                                                          _set_up_pamodel_for_simulations,
                                                          change_memprot_kcats
                                                          )
from Scripts.mcpam_generation_uniprot_id import set_up_ecolicore_mcpam


Config.BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
DATA_DIR = os.path.join('Data')  # os.path.join(os.path.split(os.getcwd())[0], 'Data')
glc_uptake_rates = list(np.linspace(1, 10, 10))


def calculate_and_plot_memprot_contribution_for_different_param(model, save_fig=False, file_number=None):
    if save_fig and file_number is None:
        raise ValueError("If save_fig=True, you must provide a file_number.")

    p_tot_list = np.linspace(0.2, 0.3, 10)
    a_available_list = np.linspace(1, 5, 5)
    data_to_plot_memprot_contribution = np.zeros((10,5))
    data_to_plot_mu = np.zeros((10,5))

    for i, p_tot in enumerate(p_tot_list):
        row_memprot_contribution = []
        row_mu = []

        for a_available in a_available_list:
            with model:
                model.p_tot = p_tot
                model.sectors.get_by_id('MembraneSector').change_available_membrane_area(a_available/100, model)
                model.optimize()
                row_memprot_contribution.append(model.sectors.get_by_id('MembraneSector').calculate_occupied_membrane(model, get_memprot_contribution = True))
                row_mu.append(model.objective.value)
        data_to_plot_memprot_contribution[i] = row_memprot_contribution
        data_to_plot_mu[i] = row_mu

    df_to_plot_memprot_contribution = pd.DataFrame(data_to_plot_memprot_contribution, columns=a_available_list, index=p_tot_list)
    df_to_plot_mu = pd.DataFrame(data_to_plot_mu, columns=a_available_list, index=p_tot_list)

    # Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Plot the first heatmap (memprot contribution to the p_tot)
    sns.heatmap(df_to_plot_memprot_contribution, annot=True, cmap='YlGnBu', cbar_kws={'label': 'Membrane Protein Contribution (%)'}, ax=axes[0])
    axes[0].set_title('Membrane Protein Contribution to Total Protein Pool')
    axes[0].set_xlabel('Available Membrane Area (units)')
    axes[0].set_ylabel('Total Protein Concentration (units)')

    # Plot the first heatmap (memprot contribution to the p_tot)
    sns.heatmap(df_to_plot_mu, annot=True, cmap='RdYlBu', cbar_kws={'label': 'Growth Rate (1/h)'}, ax=axes[1])
    axes[1].set_title('Growth Rate as a Function of Membrane Area and Protein Pool')
    axes[1].set_xlabel('Available Membrane Area (units)')
    axes[1].set_ylabel('')  # Removes the y-axis label
    axes[1].tick_params(left=False, labelleft=False)  # Hides y-axis ticks and labels

    # Show the plot
    plt.tight_layout()

    if save_fig:
        plt.savefig(f'Results/PAM_parametrizer/Simulations/2025_02_28/csc_analysis/csc_analysis_full_scale_{file_number}.png')

# ### core simulations
# mcpam_core = set_up_ecolicore_mcpam(sensitivity = False)
# mcpam_core.optimize()
# calculate_and_plot_memprot_contribution_for_different_param(mcpam_core)



### PAM simulations
#### 3.1 Build the mcPAModel
for i in range(0,7):
    file_number = i+1
    diagnostics_data_path = f'Results/PAM_parametrizer/Files/2025_02_28/pam_parametrizer_diagnostics_mciML1515_{file_number}.xlsx'
    sheet_name = 'Best_Individuals'
    pam_info_path = 'Results/PAM_parametrizer/Files/2025_02_28/proteinAllocationModel_mciML1515_EnzymaticData_multi.xlsx'

    mcpam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)
    _set_up_pamodel_for_simulations(mcpam, 'EX_glc__D_e', transl_sector_config=True)
    mcpam = create_pamodel_from_diagnostics_file(diagnostics_data_path, mcpam, sheet_name)
    calculate_and_plot_memprot_contribution_for_different_param(mcpam, save_fig=True, file_number=file_number)

    # memprot_file_path = 'Results/PAM_parametrizer/Files/2025_02_28/memprot_data.xlsx'
    # memprot_sheet_name = 'diagnostics_1'
    # mcpam = change_memprot_kcats(memprot_file_path, mcpam, memprot_sheet_name)



