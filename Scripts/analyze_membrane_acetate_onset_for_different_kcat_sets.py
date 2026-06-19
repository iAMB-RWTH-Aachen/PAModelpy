import pandas as pd
from cobra import Model
from cobra.io import read_sbml_model
from typing import Tuple, Literal
import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from src.PAModelpy.PAModel import PAModel, ActiveEnzymeSector, MembraneSector
from src.PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein, _order_enzyme_complex_id
from Scripts.mcpam_simulations_analysis import (
    run_simulation_pam_mcpam, 
    run_simulations_pam_mcpam_w_different_areas,
    run_simulations_pam_mcpam_w_different_tpcs,
    run_simulation_gem_pam_mcpam
)
from Scripts.create_pamodel_from_diagnostics_file import create_pamodel_from_diagnostics_file

def calculate_membrane_acetate_tpc_onset_difference(model, type:str = 'full scale', glc_uptake_rates:list = np.linspace(1, 10, 10)):
    # To run this function pam sensitivity has to be turned on
    if not model.sensitivity:
        model.sensitivity = True

    mus = []
    acetates = []
    membrane_sens= []
    tpc_sens = []
    delta_mus = []

    biomass_id = get_biomass_id(type=type)

    for glc in glc_uptake_rates:
        with model:
            # change glucose uptake rate
            model.change_reaction_bounds(rxn_id='EX_glc__D_e', upper_bound=-glc, lower_bound=-glc)
            # disable pyruvate formate lyase (inhibited by oxygen)
            model.change_reaction_bounds(rxn_id='PFL', upper_bound=0)
            # solve the model
            sol_pam = model.optimize()
            # save data
            mus.append(sol_pam.fluxes[biomass_id])
            acetates.append(sol_pam.fluxes['EX_ac_e'])  
            # save capacity coefficients
            capacity_coeff = model.capacity_sensitivity_coefficients
            membrane_coeff = capacity_coeff[capacity_coeff['constraint'] == 'membrane'].coefficient.to_list()
            membrane_sens.append(membrane_coeff)
            tpc_coeff = capacity_coeff[capacity_coeff['constraint'] == 'proteome'].coefficient.to_list()
            tpc_sens.append(tpc_coeff)

    mu_acetate_onset = find_onset(mus, acetates)
    mu_acetate_overflow = detect_regime_change(mus, acetates)
    mu_membrane = find_onset(mus, membrane_sens)
    mu_tpc = find_onset(mus, tpc_sens)
    delta_mu = mu_acetate_onset - mu_membrane
    
    return delta_mu, mu_membrane, mu_tpc, mu_acetate_onset, mu_acetate_overflow

def find_onset(x, y, threshold=1e-6):
    """
    Find first x value where y exceeds threshold.

    Args:
        x (array-like): List of growth rates.
        y (array-like): List of signals (acetate flux or sensitivity).
        threshold (float): Detection threshold.

    Returns:
        Onset x-value (float)
        Returns np.nan if onset is never reached.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    idx = np.where(y > threshold)[0]

    if len(idx) == 0:
        return np.nan

    return x[idx[0]]

def detect_regime_change(x, y, slope_fraction=0.3, acetate_threshold=1e-6):
    """
    Detect overflow onset as the first point where the acetate slope exceeds a fraction of the maximum observed slope

    Args:
        x (array-like): List of growth rates.
        y (array-like): List of signals (acetate flux or sensitivity).
        slope_fraction (float): Fraction of the max slope used as threshold
        acetate_threshold (float): Threshold for legit acetate value, below that is assummed to be solver noise
    Returns:
        mu (float): Growth rate at ovrflow onset
    """
    x = np.asarray(x)
    y = np.asarray(y)

    dy = np.gradient(y, x)
    threshold = slope_fraction * np.max(dy)
    idx = np.where(dy >= threshold)[0]
    if len(idx) == 0:
        return np.nan
    
    mu = x[idx[0]]

    return mu
    
def get_biomass_id(type:str):
    if type == "full scale":
        biomass_id = 'BIOMASS_Ec_iML1515_core_75p37M'
    else:
        biomass_id = 'BIOMASS_Ecoli_core_w_GAM'
    return biomass_id

def plot_onset_events(results_df):
    """
    Args:
        results_df (DataFrame):
            kcat_set
            mu_tpc
            mu_membrane
            mu_acetate
    
    Return:
        fig 
    """

    fig, ax = plt.subplots(figsize=(8, 6))

    colors = {
        "tpc": "#0072B2",
        "membrane": "#E69F00",
        "acetate_onset": "#009E73",
        "acetate_overflow": "#D55E00"
    }

    for idx, row in results_df.iterrows():
        y = idx

        # Collect the values of all columns per row
        values = [
            row['mu_tpc'],
            row['mu_membrane'],
            row['mu_acetate_onset'],
            row['mu_acetate_overflow']
        ]
        values = [v for v in values if not np.isnan(v)]

        # Plot line plot of minimum and maximum value toto connect the onset dots
        if len(values) > 1:
            ax.plot([min(values), max(values)], [y, y], linewidth=1)

        # Plot the different onsets with different symbols
        ax.scatter(row['mu_tpc'], y, marker='s', s=80,
                   label='TPC' if idx == 0 else "", color=colors["tpc"])
        ax.scatter(row['mu_membrane'], y, marker='o', s=80,
                   label='Membrane' if idx == 0 else "", color=colors["membrane"])
        ax.scatter(row['mu_acetate_onset'], y, marker='^', s=80,
                   label='Acetate onset' if idx == 0 else "", color=colors["acetate_onset"])
        ax.scatter(row['mu_acetate_overflow'], y, marker='x', s=80,
                   label='Acetate overflow' if idx == 0 else "", color=colors["acetate_overflow"])

    ax.set_yticks(range(len(results_df)))
    ax.set_yticklabels(results_df['kcat_set'])

    ax.set_xlabel('Growth rate μ [h$^{-1}$]')
    ax.set_ylabel('kcat set')
    ax.legend()

    return fig

if __name__ == '__main__':
    pam_info_file = 'Results/PAM_parametrizer/Diagnostics_files/2026_06_10/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
    model_path = 'Models/iML1515.xml'

    mcpam = set_up_pam(pam_info_file=pam_info_file, 
                       model=model_path,
                       sensitivity=True, 
                       membrane_sector=True,
                       separate_memprot_from_tpc=True,
                       total_protein=0.241,
                       usable_area_fraction=0.6174,
                       enable_unused_membrane_sector=True
                       )
    results = []
    for file_nr in range(1, 10):
        with mcpam:
            mcpam = create_pamodel_from_diagnostics_file(file_path=f'Results/PAM_parametrizer/Diagnostics_files/2026_06_10/pam_parametrizer_diagnostics_mciML1515_{file_nr}.xlsx',
                                                    model=mcpam,
                                                    sheet_name='Best_Individuals')

            delta_mu, mu_membrane, mu_tpc, mu_acetate_onset, mu_acetate_overflow = calculate_membrane_acetate_tpc_onset_difference(mcpam)
            kcat_set = f"2026_06_10_{file_nr}"
            results.append([kcat_set, mu_membrane, mu_tpc, mu_acetate_onset, mu_acetate_overflow])
    
    results_df = pd.DataFrame(results, columns=['kcat_set', 'mu_membrane', 'mu_tpc', 'mu_acetate_onset', 'mu_acetate_overflow'])
    fig = plot_onset_events(results_df)

    os.makedirs("Figures", exist_ok=True)
    fig.savefig(
        f"Figures/acetate_onset_analysis_2026_06_10.png",
        dpi=250,
        bbox_inches='tight'
    )


    