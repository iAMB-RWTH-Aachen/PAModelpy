from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import seaborn as sns
import re

import pandas as pd
import os
import sys
import numpy as np

sys.path.append('C:\\Users\\claud\\Documents\\iamb-student-folders\\iamb-folder-template\\mcPAM_package')

from src.PAModelpy.configuration import Config
from src.PAModelpy.MembraneSector import UnusedMembraneSector

if os.path.split(os.getcwd())[1] == 'Figures':
    os.chdir(os.path.split(os.getcwd())[0])
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
from Scripts.create_pamodel_from_diagnostics_file import (create_pamodel_from_diagnostics_file,
                                                          _set_up_pamodel_for_simulations,
                                                          )


Config.BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
DATA_DIR = os.path.join('Data')  # os.path.join(os.path.split(os.getcwd())[0], 'Data')

### 1 Useful functions
def calculate_sensitivities(pamodel, glc_uptake_rates:list = None):
    glc_uptake_rates = np.linspace(1, 10, 10) if glc_uptake_rates is None else glc_uptake_rates
    Ccsc = []
    Cesc = []
    y_axis = []
    fluxes = []

    # disable pyruvate formate lyase (inhibited by oxygen)

    for glc in glc_uptake_rates:
        print('glucose uptake rate ', glc, ' mmol/gcdw/h')
        with pamodel:
            # change glucose uptake rate
            pamodel.change_reaction_bounds(rxn_id='EX_glc__D_e', upper_bound=-glc, lower_bound=-glc)
            pamodel.change_reaction_bounds(rxn_id='PFL', upper_bound=0)
            # solve the model
            sol_pam = pamodel.optimize()
            fluxes.append(sol_pam.fluxes)
            if pamodel.solver.status == 'optimal': y_axis += [glc] 
            # save data
            Ccsc_new = list()

            if pamodel.solver.status == 'optimal':
                capacity_coeff = pamodel.capacity_sensitivity_coefficients
                # for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']: # for PAM without membrane
                for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector', 'membrane']:
                    Ccsc_new += capacity_coeff[capacity_coeff['constraint'] == csc].coefficient.to_list()

                Ccsc += [Ccsc_new]

                enzyme_coeff = pamodel.enzyme_sensitivity_coefficients
                Cesc += [enzyme_coeff.coefficient.to_list()]

                print('Sum of capacity sensitivity coefficients: \t \t \t \t \t \t', round(sum(Ccsc_new), 6))
                print('Sum of enzyme sensitivity coefficients: \t \t \t \t \t \t', round(sum(Cesc[-1]), 6), '\n')

    return {'Ccsc': Ccsc, 'Cesc': Cesc, 'y_axis': y_axis, 'fluxes': fluxes, 'capacity coefficients': capacity_coeff,
            'enzyme coefficients': enzyme_coeff}


# %%
def parse_x_axis_heatmap(capacity_coeff, enzyme_coeff):
    x_axis_csc = []

    # for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']: # for PAM without membrane
    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector', 'membrane']:
        if csc == 'flux_ub' or csc == 'flux_lb':
            x_axis_csc += [coef + '_' + csc for coef in
                           capacity_coeff[capacity_coeff['constraint'] == csc].rxn_id.to_list()]
        else:
            x_axis_csc += [coef + '_' + csc for coef in capacity_coeff[
                capacity_coeff['constraint'] == csc].enzyme_id.to_list()]

    x_axis_esc = enzyme_coeff.enzyme_id.to_list()
    return x_axis_csc, x_axis_esc


# %%
def make_heatmap_subfigure(results, csc_matrix, esc_matrix, x_csc, x_esc, yaxis, fig, grdspc,
                           ylabels=True, xlabels=False, cbar=True, title=None,
                           fontsize=16, vmin=-1.5, vmax=1.5,
                           annotate=None, phenotype_data=None, cmap=None):

    if cmap is None:
        colors_neg = plt.cm.Blues(np.linspace(1, 0.3, 128))
        colors_pos = plt.cm.OrRd(np.linspace(0.1, 1, 128))
        colors_zero = np.array([[1, 1, 1, 1]])

        colors = np.vstack((colors_neg, colors_zero, colors_pos))
        combined_cmap = mcolors.ListedColormap(colors)

        bounds = np.linspace(vmin, vmax, len(colors))
        norm = mcolors.BoundaryNorm(bounds, combined_cmap.N)
    else:
        combined_cmap = cmap
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    if cbar:
        gs = gridspec.GridSpecFromSubplotSpec(
            3, 2,
            width_ratios=[len(yaxis), 1],
            height_ratios=[1, len(x_csc), len(x_esc)],
            hspace=0,
            subplot_spec=grdspc
        )
    else:
        gs = gridspec.GridSpecFromSubplotSpec(
            3, 1,
            width_ratios=[len(yaxis)],
            height_ratios=[1, len(x_csc), len(x_esc)],
            hspace=0,
            subplot_spec=grdspc
        )

    esc_ax = fig.add_subplot(gs[2, 0])
    acetate_ax = fig.add_subplot(gs[0, 0], sharex=esc_ax)
    csc_ax = fig.add_subplot(gs[1, 0], sharex=esc_ax)

    if cbar:
        cbar_ax = fig.add_subplot(gs[1:, 1])

    if annotate is not None:
        acetate_ax.annotate(
            annotate,
            xy=(2, 1),
            xycoords="data",
            xytext=(-0.05, 1.5),
            textcoords="axes fraction",
            va="top",
            ha="left",
            fontsize=fontsize * 1.5,
            weight="bold"
        )

    glc_fluxes = np.array([-sim.EX_glc__D_e for sim in results["fluxes"]])
    acetate_fluxes = np.array([sim.EX_ac_e for sim in results["fluxes"]])

    x_start = glc_fluxes.min()
    x_end = glc_fluxes.max()

    glc_onset = None
    for g, ac in zip(glc_fluxes, acetate_fluxes):
        if ac > 0.01:
            glc_onset = g
            break

    acetate_ax.spines["top"].set_visible(False)
    acetate_ax.spines["right"].set_visible(False)
    
    # ------------------------------
    # Respiration / overflow arrows
    # ------------------------------

    y_arrow = 11

    if glc_onset is None:
        # Entire range is respiration
        acetate_ax.annotate(
            "",
            xy=(x_end, y_arrow),
            xytext=(x_start, y_arrow),
            arrowprops=dict(arrowstyle="<->", color="purple", lw=2),
        )

        acetate_ax.text(
            (x_start + x_end) / 2,
            y_arrow + 0.8,
            "Respiration",
            color="purple",
            ha="center",
            va="bottom",
            fontsize=fontsize,
        )

    else:
        # Respiration
        acetate_ax.annotate(
            "",
            xy=(glc_onset, y_arrow),
            xytext=(x_start, y_arrow),
            arrowprops=dict(arrowstyle="<->", color="purple", lw=2),
        )

        acetate_ax.text(
            (x_start + glc_onset) / 2,
            y_arrow + 0.8,
            "Respiration",
            color="purple",
            ha="center",
            va="bottom",
            fontsize=fontsize,
        )

        # Overflow
        acetate_ax.annotate(
            "",
            xy=(x_end, y_arrow),
            xytext=(glc_onset, y_arrow),
            arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        )

        acetate_ax.text(
            (glc_onset + x_end) / 2,
            y_arrow + 0.8,
            "Overflow",
            color="black",
            ha="center",
            va="bottom",
            fontsize=fontsize,
        )

    # ------------------------------
    # Acetate production curve
    # ------------------------------

    acetate_ax.plot(
        glc_fluxes,
        acetate_fluxes,
        linewidth=4,
        color="darkblue",
    )

    if phenotype_data is not None:
        acetate_ax.scatter(
            phenotype_data["EX_glc__D_e"],
            phenotype_data["EX_ac_e"],
            color="purple",
            marker="o",
            s=40,
            linewidths=1.3,
            zorder=3,
            label="Data",
        )

    acetate_ax.set_xlim(
        glc_fluxes.min() - 0.25,
        glc_fluxes.max() + 0.25,
    )

    acetate_ax.set_ylim(-0.5, max(15, acetate_fluxes.max() + 1))

    acetate_ax.tick_params(axis="y", labelsize=fontsize)
    acetate_ax.xaxis.set_visible(False)

    if ylabels:
        acetate_ax.set_ylabel(
            r"Acetate"
            "\n"
            r"[$mmol_{ac}/g_{CDW}/h$]",
            fontsize=fontsize,
            rotation=0,
        )

    if title is not None:
        acetate_ax.set_title(title, fontsize=fontsize * 1.5)

    # ------------------------------
    # Heatmap extent
    # ------------------------------

    # Compute heatmap cell boundaries from the glucose uptake values.
    if len(yaxis) > 1:
        dx = np.diff(yaxis).mean()
    else:
        dx = 1.0

    xmin = yaxis[0] - dx / 2
    xmax = yaxis[-1] + dx / 2

    # ------------------------------
    # CSC heatmap
    # ------------------------------

    im_csc = csc_ax.imshow(
        csc_matrix,
        aspect="auto",
        cmap=combined_cmap,
        norm=norm,
        extent=[xmin, xmax, len(x_csc) - 0.5, -0.5],
    )

    csc_ax.set_yticks(np.arange(len(x_csc)))
    csc_ax.set_yticklabels(x_csc, fontsize=fontsize)

    csc_ax.xaxis.set_visible(False)

    if ylabels:
        csc_ax.set_ylabel("CSC", fontsize=fontsize * 1.25)

    csc_ax.spines["bottom"].set_linewidth(10)
    csc_ax.spines["bottom"].set_color("black")

    # ------------------------------
    # ESC heatmap
    # ------------------------------

    im_esc = esc_ax.imshow(
        esc_matrix,
        aspect="auto",
        cmap=combined_cmap,
        norm=norm,
        extent=[xmin, xmax, len(x_esc) - 0.5, -0.5],
    )

    esc_ax.set_yticks(np.arange(len(x_esc)))
    esc_ax.set_yticklabels(x_esc, fontsize=fontsize)

    # IMPORTANT:
    # Use the actual glucose uptake values instead of column numbers.
    esc_ax.set_xticks(yaxis)
    esc_ax.set_xticklabels(
        [str(v) for v in yaxis],
        rotation=45,
        ha="right",
        fontsize=fontsize,
    )

    if ylabels:
        esc_ax.set_ylabel("ESC", fontsize=fontsize * 1.25)

    if xlabels:
        esc_ax.set_xlabel(
            r"Glucose uptake rate [$mmol_{glc}/g_{CDW}/h$]",
            fontsize=fontsize * 1.25,
        )

    # ------------------------------
    # Colorbar
    # ------------------------------

    if cbar:
        cbar_ax.xaxis.set_visible(False)

        make_scaled_colorbar(
            ax=cbar_ax,
            fig=fig,
            cmap=combined_cmap,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            fontsize=fontsize * 1.25,
        )

    fig.set_figwidth(24)
    fig.set_figheight(7)

    fig.align_labels()

    return fig

def make_scaled_colorbar(ax, fig, cmap, norm, vmin, vmax,
                         fontsize=16, cbarlabel='Sensitivity Coefficient'):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax, cax=ax, shrink=0.2, fraction=1)

    # Adjust the tick intervals
    tick_locations = np.linspace(vmin, vmax, num=5)  # Adjust num to the desired number of ticks
    cbar.set_ticks(tick_locations)
    cbar.set_ticklabels([f"{tick:.1f}" for tick in tick_locations])  # Optional: customize tick labels

    # Setting the fontsize of the colorbar
    cbar.set_label(cbarlabel, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.yaxis.get_offset_text().set(size=fontsize)


# %%
# adjust labels for better readibility
def adjust_heatmap_labels(labels):
    new_labels = labels.copy()
    for i, label in enumerate(labels):
        if 'EX_glc__D_e' in label or label[:-3] == 'EX_glc__D_e':
            if label[-1] == 'B':
                new_labels[i] = 'EX_glc_' + label[-2:]
            else:
                new_labels[i] = 'Glucose uptake lb'

        if label == 'ATPM_flux_lb':
            new_labels[i] = "ATP maintanance lb"

        if label == 'TotalProteinConstraint_proteome':
            new_labels[i] = 'Protein pool'

        if label == 'MembraneSector_membrane':
            new_labels[i] = 'Membrane sector'

        if label[-1].isdigit() and len(label) > 3:  # all enzyme ids start with a digit

            rxn_ids = mcpam.get_reactions_with_enzyme_id(label)
            id = rxn_ids[0].split('_')
            rxn_name = id[1]
            if rxn_name in new_labels:
                new_labels[i] = rxn_name + f'*'
            else:
                new_labels[i] = '\n'.join([part for part in rxn_name.split(' ')])

        # if label[-1].isdigit() and len(label) == 3:
        #     rxn_ids = mcpam_core.get_reactions_with_enzyme_id(label)
        #     rxn_name = mcpam_core.reactions.get_by_id(rxn_ids[-1]).name.split('(')[0]
        #     new_labels[i] = '\n'.join([part for part in rxn_name.split(' ')])

            # if len(rxn_ids)>2:
            #     new_labels[i] = pamodel.reactions.get_by_id(rxn_ids[-1]).name.split('(')[0]
            # else:

            #     new_labels[i] = ',\n'.join([pamodel.reactions.get_by_id(rxn_ids[-1]).name.split('(')[0]for rxn in rxn_ids])
    return new_labels


# %%
def find_nonzero_sensitivities(Cv, x_axis):
    indices = []
    for row in Cv:
        for index, coeff in enumerate(row):
            if abs(coeff) > 0 and index not in indices:
                indices.append(index)

    coeff_nonzero = []
    for row in Cv:
        coeff_nonzero.append([coeff for i, coeff in enumerate(row) if i in indices])
    x_coeff_nonzero = [coeff for i, coeff in enumerate(x_axis) if i in indices]

    return coeff_nonzero, x_coeff_nonzero


# %%
def find_top5_sensitivities(Cv, x_axis, yaxis, threshold=0.01):
    # top 5 enzymes per simulation
    Cv_df = pd.DataFrame(Cv, columns=x_axis, index=yaxis)
    largest = list()
    for i, row in Cv_df.iterrows():
        top5 = abs(row).nlargest()
        if top5.iloc[0]:
            largest += [index for index, value in top5.items() if abs(value) > threshold]
    # remove duplicates
    largest_list = list(set(largest))

    # extract non duplicate top5 enzymes
    top5_df = Cv_df[largest_list].T.drop_duplicates().sort_index()
    largest_list = top5_df.index.values

    top5_matrix = [list(row) for i, row in top5_df.iterrows()]
    return top5_matrix, largest_list

##############################################################################
# BUILD MODEL
##############################################################################

pam_info_path = 'Results/PAM_parametrizer/Diagnostics_files/2026_06_10/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
model_path = 'Models/iML1515.xml'
diagnostics_path = 'Results/PAM_parametrizer/Diagnostics_files/2026_06_10/pam_parametrizer_diagnostics_mciML1515_4.xlsx'

mcpam = set_up_pam(pam_info_file=pam_info_path, 
                    model=model_path,
                    sensitivity=True, 
                    membrane_sector=True,
                    separate_memprot_from_tpc=True,
                    total_protein=0.241,
                    usable_area_fraction=0.6174,
                    unused_membrane_sector=UnusedMembraneSector()
                    )

area = mcpam.sectors.get_by_id('MembraneSector').usable_area_fraction

mcpam = create_pamodel_from_diagnostics_file(file_path=diagnostics_path,
                                                 model=mcpam,
                                                 sheet_name='Best_Individuals')
mcpam.change_sector_parameters(mcpam.unused_enzymes, slope=0.012, lin_rxn_id='EX_glc__D_e', print_change=True)
mcpam.sectors.get_by_id('MembraneSector').unused_membrane_sector.change_unused_membrane_parameters(unused_membrane_fraction=0.12)
suffix = diagnostics_path.split('_')[-1].split('.')[0]
##############################################################################
# PHENOTYPE DATA
##############################################################################

pt_data = pd.read_excel(
    os.path.join(DATA_DIR, 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls'),
    sheet_name='Yields'
)
pt_data['EX_glc__D_e'] = -pt_data['EX_glc__D_e']

##############################################################################
# PARAMETERS
##############################################################################

fontsize = 26
width = 40
height = 14
glc_uptake_rates = np.arange(1, 10, 0.5)

##############################################################################
# MAIN LOOP: ONE FIGURE PER AREA
##############################################################################

print(f"\nRunning mcPAM with {area*100:.1f}% membrane area\n")

results = calculate_sensitivities(mcpam, glc_uptake_rates=glc_uptake_rates)

# Parse axes
x_axis_csc, x_axis_esc = parse_x_axis_heatmap(
    results['capacity coefficients'],
    results['enzyme coefficients']
)

# CSC → nonzero
csc_nonzero, x_csc_nonzero = find_nonzero_sensitivities(
    results['Ccsc'], x_axis_csc
)
csc_nonzero_t = np.transpose(np.array(csc_nonzero))

# ESC → top5
esc_top5, x_esc_top5 = find_top5_sensitivities(
    results['Cesc'],
    x_axis_esc,
    glc_uptake_rates
)
esc_top5_t = np.transpose(np.array(esc_top5))

# Adjust labels
x_csc_labels = adjust_heatmap_labels(x_csc_nonzero)
x_esc_labels = adjust_heatmap_labels(x_esc_top5)

##########################################################################
# CREATE FIGURE (NO OUTER GRIDSPEC)
##########################################################################
fig = plt.figure(layout='constrained')
gs = fig.add_gridspec(1, 1)
gs_pam = gs[0]

make_heatmap_subfigure(
    results=results,
    csc_matrix=csc_nonzero_t,
    esc_matrix=esc_top5,
    x_csc=x_csc_labels,
    x_esc=x_esc_labels,
    yaxis=glc_uptake_rates,
    fig=fig,
    grdspc=gs_pam,
    ylabels=True,
    xlabels=True,
    phenotype_data=pt_data,
    fontsize=fontsize,
    title=f"mcPAM – {area*100:.1f}% membrane area"
)

fig.set_figwidth(width)
fig.set_figheight(height)
os.makedirs("Figures", exist_ok=True)
fig.savefig(
    f"Figures/sENZ_2026_06_10/mcPAM_sensitivities_diagnostics_{suffix}_unused_fra_0_12_ue_slope_12mg.png",
    dpi=250,
    bbox_inches='tight'
)

plt.close(fig)
