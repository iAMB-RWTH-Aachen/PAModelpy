import cobra
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# load PAMpy modules
if os.path.split(os.getcwd())[1] == 'Scripts':
    os.chdir('..')

from src.PAModelpy import PAModel, ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.PAModelpy.PAMValidator import PAMValidator
from Scripts.pam_generation_uniprot_id import set_up_ecolicore_pam, set_up_ecolicore_mcpam


mcpam_core = set_up_ecolicore_mcpam(sensitivity=False)

#### Fluxes simulations for different glc uptake rates
# load phenotype data from excel file
pt_data = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes','Ecoli_phenotypes_py_rev.xls'),
                        sheet_name='Yields', index_col=None)

# extract reaction specific data
rxn_to_pt = {}
rxn_transform = {
    'EX_ac_e': 'EX_ac_e',
    'EX_co2_e': 'EX_co2_e',
    'EX_o2_e': 'EX_o2_e',
    'BIOMASS_Ecoli_core_w_GAM':'BIOMASS_Ec_iML1515_core_75p37M'
    # 'BIOMASS_Ecoli_core_w_GAM':'BIOMASS_Ecoli_core_w_GAM'
}
for rxn_id, pt_id in rxn_transform.items():
    rxn_to_pt[rxn_id] = pt_data[['EX_glc__D_e', pt_id]].dropna().rename(columns={pt_id: rxn_id})

with mcpam_core:
    # change glucose uptake rate
    mcpam_core.reactions.EX_glc__D_e.lower_bound = -6.0
    # solve the model
    sol_pam = mcpam_core.optimize()
    # print(pamodel.summary())
    # with pd.option_context('display.max_rows', None):
    #     print(sol_pam.fluxes)
mcpam_core.optimize()

glc_uptake_rates = np.linspace(0.5, 11.5, 20)
fluxes = []
concentrations = [0]
for glc in glc_uptake_rates:
    with mcpam_core:
        # change glucose uptake rate
        mcpam_core.reactions.EX_glc__D_e.lower_bound = -glc
        # disable pyruvate formate lyase (inhibited by oxygen)
        mcpam_core.reactions.PFL.upper_bound = 0
        # solve the model
        sol_pam = mcpam_core.optimize()
        # save data
        fluxes.append(sol_pam.fluxes)  # flux distributions
        concentration = 0
        for enz_var in mcpam_core.enzyme_variables:
            concentration += enz_var.concentration
        concentrations.append(concentration)

# plot flux changes with glucose uptake
rxn_id = ['EX_ac_e', 'EX_co2_e', 'EX_o2_e', 'BIOMASS_Ecoli_core_w_GAM']
fig, axs = plt.subplots(2, 2, dpi=90)
for r, ax in zip(rxn_id, axs.flatten()):
    # plot data
    if r in rxn_to_pt.keys():
        ax.scatter(abs(rxn_to_pt[r]['EX_glc__D_e']), abs(rxn_to_pt[r][r]),
                   color='firebrick', marker='o', s=30, linewidths=1.3,
                   facecolors=None, zorder=0,
                   label='Data')

    # plot simulation
    ax.plot(glc_uptake_rates, [abs(f[r]) for f in fluxes],
            label='Simulation', linewidth=2.5,
            zorder=5)

    # options
    ax.set_xlabel('glc uptake rate [mmol/gDW/h]')
    ax.set_ylabel('flux [mmol/gDW/h]')
    ax.set_title(r)
    # set grid
    ax.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)
    # show legend
    ax.legend(fontsize=8, edgecolor='white', facecolor='white', framealpha=1)

# Add parameter box
param_text = (
    f"mcPAM_core model \n"
    f"Total protein: 0.16995 g/g DW \n"
    f"Max membrane area = 28.5%"
)
# Position the text box in the upper left corner of each subplot
fig.text(0.6, 2.7, param_text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, boxstyle="round,pad=0.3"))

plt.tight_layout()
plt.show()
