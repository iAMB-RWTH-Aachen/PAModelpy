from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors

import pandas as pd
import os
import sys
import numpy as np

sys.path.append('C:\\Users\\claud\\Documents\\iamb-student-folders\\iamb-folder-template\\mcPAM_package')

from src.PAModelpy.configuration import Config
if os.path.split(os.getcwd())[1] == 'Figures':
    os.chdir(os.path.split(os.getcwd())[0])
from Scripts.mcpam_generation_uniprot_id import set_up_ecolicore_pam, set_up_ecoli_mcpam
from Scripts.mcpam_simulations_analysis import run_pam_mcpam_core_with_optimized_kcats

Config.BIOMASS_REACTION = 'BIOMASS_Ec_iML1515_core_75p37M'
DATA_DIR = os.path.join('Data')  # os.path.join(os.path.split(os.getcwd())[0], 'Data')
glc_uptake_rates = list(np.linspace(1, 10, 10))

### 1 Useful functions
def calculate_sensitivities(pamodel):
    glc_uptake_rates = np.linspace(1, 10, 10)
    Ccsc = []
    Cesc = []
    y_axis = []
    fluxes = []

    # disable pyruvate formate lyase (inhibited by oxygen)
    pamodel.change_reaction_bounds(rxn_id='PFL', upper_bound=0)

    for glc in glc_uptake_rates:
        print('glucose uptake rate ', glc, ' mmol/gcdw/h')
        with pamodel:
            # change glucose uptake rate
            pamodel.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                           lower_bound=-glc, upper_bound=-glc)
            # solve the model
            # pamodel.objective = 'EX_ac_e'
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
                           ylabels=True, xlabels=False, cbar=True, title=None, fontsize=5,
                           vmin=-1.5, vmax=1.5, annotate=None, phenotype_data=None, cmap=None
                           # cmap = plt.cm.get_cmap('viridis')
                           ):
    # fig = plt.figure()
    # adjust labels for better readibility
    x_csc = adjust_heatmap_labels(x_csc_nonzero_pam)
    x_esc = adjust_heatmap_labels(x_esc_top5_pam)

    if cmap is None:
        # Create separate colormaps for positive and negative values and a color for zero
        colors_neg = plt.cm.Blues(np.linspace(1, 0.3, 128))
        colors_pos = plt.cm.OrRd(np.linspace(0.1, 1, 128))  # plt.cm.Reds(np.linspace(0, 0.5, 128))
        colors_zero = np.array([[1, 1, 1, 1]])  # gray for zero

        # Combine them into a single colormap
        colors = np.vstack((colors_neg, colors_zero, colors_pos))
        combined_cmap = mcolors.ListedColormap(colors, name='custom_cmap')

        # Create a norm that handles the zero color properly
        bounds = np.linspace(vmin, vmax, len(colors))
        norm = mcolors.BoundaryNorm(bounds, combined_cmap.N)

    if cbar:
        gs = gridspec.GridSpecFromSubplotSpec(3, 2, width_ratios=[len(yaxis), 1],
                                              height_ratios=[1, len(x_csc), len(x_esc)], hspace=0,
                                              subplot_spec=grdspc)
    else:
        gs = gridspec.GridSpecFromSubplotSpec(3, 1, width_ratios=[len(yaxis)],
                                              height_ratios=[1, len(x_csc), len(x_esc)], hspace=0,
                                              subplot_spec=grdspc)

    acetate_ax = fig.add_subplot(gs[0, 0])  # acetate production
    csc_ax = fig.add_subplot(gs[1, 0])  # CSC heatmap
    esc_ax = fig.add_subplot(gs[2, 0], sharex=csc_ax)  # ESC heatmap

    if cbar:
        cbar_ax = fig.add_subplot(gs[1:, 1])  # colorbar

    # add annotation for subfigure (A or B)
    if annotate is not None:
        acetate_ax.annotate(annotate, xy=(2, 1), xycoords='data',
                            xytext=(-0.05, 1.5), textcoords='axes fraction',
                            va='top', ha='left', fontsize=fontsize * 1.5, weight='bold')

    glc_fluxes = [-sim.EX_glc__D_e for sim in results['fluxes']]

    # add arrow indicating growth regime
    # 0. remove the box to improve readability of the text
    acetate_ax.spines['top'].set_visible(False)
    acetate_ax.spines['right'].set_visible(False)

    # 1. Find the start of the overflow regime (which is when acetate is being produced)
    for i, ac in enumerate([sim.EX_ac_e for sim in results['fluxes']]):
        if ac > 0.001:
            glc_onset = glc_fluxes[i] - 1
            break
    # 2. determine the dx covered by respiration
    dx_respiration = glc_onset - glc_fluxes[0]
    # 3 create respiration arrow
    # forward arrow
    acetate_ax.arrow(
        glc_fluxes[0], 11, dx_respiration, 0,
        linewidth=2, color='purple', label='Respiration', length_includes_head=True, head_width=3, head_length=0.5
    )
    # reverse arrow
    acetate_ax.arrow(
        glc_onset, 11, -dx_respiration, 0, head_starts_at_zero=True,
        linewidth=2, color='purple', label='Respiration', length_includes_head=True, head_width=3, head_length=0.5
    )
    # annotate
    acetate_ax.annotate('Respiration',
                        xy=(dx_respiration / 3, 15),
                        xytext=(10, -10), fontsize=fontsize,
                        textcoords='offset points', color='purple')
    # remove the box
    acetate_ax.spines['top'].set_visible(False)
    acetate_ax.spines['right'].set_visible(False)

    # 4. create overflow arrow
    # forward arrow
    acetate_ax.arrow(
        glc_onset, 11, 10 - glc_onset, 0,
        linewidth=2, color='black', label='Overflow', length_includes_head=True, head_width=3, head_length=0.5
    )
    # reverse arrow
    acetate_ax.arrow(
        10, 11, -(10 - glc_onset), 0,
        linewidth=2, color='black', label='Overflow', length_includes_head=True, head_width=3, head_length=0.5
    )
    # annotate
    acetate_ax.annotate('Overflow', fontsize=fontsize,
                        xy=((10 - glc_onset - 2) / 2 + glc_onset, 15),
                        xytext=(10, -10),
                        textcoords='offset points', color='black')

    # acetate graph
    acetate_ax.plot([-sim.EX_glc__D_e for sim in results['fluxes']], [sim.EX_ac_e for sim in results['fluxes']],
                    linewidth=4, color='darkblue')
    acetate_ax.tick_params(axis='y', labelsize=fontsize)
    acetate_ax.set_xlim([0, 10.5])
    acetate_ax.set_ylim([-0.5, 15])
    acetate_ax.xaxis.set_visible(False)
    if ylabels:
        acetate_ax.set_ylabel(r'Acetate' '\n' '[$mmol_{ac}/g_{CDW}/h$]', fontsize=25, rotation=0, position=(1,-5))

    # add phenotype data if this is given
    if phenotype_data is not None:
        acetate_ax.scatter(phenotype_data['EX_glc__D_e'], phenotype_data['EX_ac_e'],
                           color='purple', marker='o', s=40, linewidths=1.3,
                           facecolors=None, zorder=0,
                           label='Data')

    if title is not None: acetate_ax.set_title(title, fontsize=fontsize * 1.5)

    # CAC heatmap
    im_csc = csc_ax.imshow(csc_matrix, aspect="auto", cmap=combined_cmap, norm=norm)
    csc_ax.set_yticks(np.arange(len(x_csc)), labels=x_csc, fontsize=fontsize)
    csc_ax.xaxis.set_visible(False)
    if ylabels:
        csc_ax.set_ylabel('CSC', fontsize=fontsize * 1.25)

    # Make line between CSC and ESC data more clear
    axis = 'bottom'
    csc_ax.spines[axis].set_linewidth(10)
    csc_ax.spines[axis].set_color("black")
    csc_ax.spines[axis].set_zorder(0)

    # ESC heatmap
    im_esc = esc_ax.imshow(esc_matrix, aspect="auto", cmap=combined_cmap, norm=norm)
    esc_ax.set_yticks(np.arange(len(x_esc)), labels=x_esc, fontsize=fontsize)
    esc_ax.set_xticks(np.arange(len(yaxis)), labels=yaxis, fontsize=fontsize, rotation=45, ha='right')
    if ylabels:
        esc_ax.set_ylabel('ESC', fontsize=fontsize * 1.25)
    if xlabels:
        esc_ax.set_xlabel('Glucose uptake rate [$mmol_{glc}/g_{CDW}/h$]', fontsize=fontsize * 1.25)

    # colorbar
    if cbar:
        cbar_ax.xaxis.set_visible(False)
        make_scaled_colorbar(ax=cbar_ax, fig=fig, cmap=combined_cmap, norm=norm,
                             vmin=vmin, vmax=vmax, fontsize=fontsize * 1.25)
    fig.set_figwidth(24)
    fig.set_figheight(7)
    fig.align_labels()
    return fig


# %%

def make_scaled_colorbar(ax, fig, cmap, norm, vmin, vmax,
                         fontsize=16, cbarlabel='Sensitivity Coefficient'):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax, cax=ax, shrink=1, fraction=1)

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
                new_labels[i] = 'EX_glc_lb'
        if label == 'TotalProteinConstraint_proteome':
            new_labels[i] = 'Protein pool'

        if label == 'MembraneSector_membrane':
            new_labels[i] = 'Membrane sector'

        if label[-1].isdigit() and len(label) > 3:  # all enzyme ids start with a digit
            rxn_ids = mcpam.get_reactions_with_enzyme_id(label)

            id = rxn_ids[0].split('_')
            rxn_name = id[1]
            new_labels[i] = '\n'.join([part for part in rxn_name.split(' ')])
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

### PAM simulations
#### 3.1 Build PAModel

pam, mcpam = run_pam_mcpam_core_with_optimized_kcats(type="full scale")

#### 3.2 Run simulations for glucose uptake of 0-10 mmol/gcdw/h
# %%
results_pam = calculate_sensitivities(mcpam)
x_axis_csc_pam, x_axis_esc_pam = parse_x_axis_heatmap(results_pam['capacity coefficients'],
                                                      results_pam['enzyme coefficients'])
# %%
# get nonzero sensitivities
csc_nonzero_pam, x_csc_nonzero_pam = find_nonzero_sensitivities(results_pam['Ccsc'], x_axis=x_axis_csc_pam)
esc_nonzero_pam, x_esc_nonzero_pam = find_nonzero_sensitivities(results_pam['Cesc'], x_axis=x_axis_esc_pam)
csc_nonzero_pam_t = np.transpose(np.array(csc_nonzero_pam))
esc_nonzero_pam_t = np.transpose(np.array(esc_nonzero_pam))

# get top5 nonzero sensitivities
csc_top5_pam, x_csc_top5_pam = find_top5_sensitivities(results_pam['Ccsc'], x_axis=x_axis_csc_pam,
                                                       yaxis=glc_uptake_rates)
esc_top5_pam, x_esc_top5_pam = find_top5_sensitivities(results_pam['Cesc'], x_axis=x_axis_esc_pam,
                                                       yaxis=glc_uptake_rates)
csc_top5_pam_t = np.transpose(np.array(csc_top5_pam))
esc_top5_pam_t = np.transpose(np.array(esc_top5_pam))

### 4 Create plot

#### 4.1 Load phenotypic data

# load phenotype data from excel file
pt_data = pd.read_excel(os.path.join(DATA_DIR, 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls'), sheet_name='Yields',
                        index_col=None)
pt_data['EX_glc__D_e'] = -pt_data['EX_glc__D_e']

# create 2 plots: supplements and main text
fontsize = 28
width = 50
height = 10
# select colormap
cmap = None  # plt.cm.get_cmap('magma')

# gridspec inside gridspec
fig = plt.figure(layout='constrained')

gs0 = gridspec.GridSpec(1, 1, figure=fig)
gs_pam = gs0[0]

fig_pam = make_heatmap_subfigure(results=results_pam, csc_matrix=csc_nonzero_pam_t, esc_matrix=esc_top5_pam,
                                 ylabels=True, xlabels=True, x_csc=x_csc_nonzero_pam, x_esc=x_esc_top5_pam,
                                 yaxis=glc_uptake_rates, fig=fig, grdspc=gs_pam,
                                 phenotype_data=pt_data, fontsize=fontsize, cmap=cmap)
fig_pam.subplots_adjust(left=0.3)
# set common x axis title
# fig_pam.xlabel('Glucose uptake rate [$mmol_{glc}/g_{CDW}/h$]', fontsize = fontsize*1.25)

# # add image
#
# ax_fig = fig.add_subplot(gs_figure)
# ax_fig.imshow(sensitivities_mapped)
# ax_fig.annotate('B', xy=(2, 1), xycoords='data',
#                 xytext=(-0.05, 1.30), textcoords='axes fraction',
#                 va='top', ha='left', fontsize=fontsize * 1.5, weight='bold')
# ax_fig.axis('off')
# ax_fig.set_xticks([])
# ax_fig.set_yticks([])

plt.plasma()
fig.subplots_adjust(left=0.3)
fig.set_figwidth(width)
fig.set_figheight(height)
fig.align_labels()

plt.show()

# fig.savefig('Figures/Figure_sensitivities_pam_uniprotid_(2).png', dpi=200, bbox_inches='tight')
