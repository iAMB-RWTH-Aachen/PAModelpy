import numpy as np
import pandas as pd
from cobra.io import load_json_model
import plotly.express
import math
from PIL import Image
import os
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors

from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector
from src.PAModelpy.MembraneSector import MembraneSector
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.configuration import Config


def build_toy_model(sensitivity:bool=True):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'R11'
    config.GLUCOSE_EXCHANGE_RXNID = 'R2'
    config.CO2_EXHANGE_RXNID = 'R12'
    config.ACETATE_EXCRETION_RXNID = 'R7'

    nmbr_reactions = 12

    # Building Active Enzyme Sector
    kcat_fwd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 10e-2, 10e-2]
    kcat_rev = [kcat for kcat in kcat_fwd]
    gpr_string = [['gene1'], ['gene2'],['gene3'],['gene4'],['gene5'],['gene6'],['gene7'],['gene8']]
    pra_string = ['E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9']

    rxn2protein = {}
    protein2gene = {}
    for i in range(nmbr_reactions-4): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+2}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass
        rxn2protein = {**rxn2protein, **{rxn_id: {f'E{i+2}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6),
                                                             'molmass': 1e6, 'genes': gpr_string[i],
                                                             'protein_reaction_association': pra_string[i]}}}}
        protein2gene = {**protein2gene, **{f'E{i+2}': gpr_string[i]}}

    active_enzyme = ActiveEnzymeSector(rxn2protein = rxn2protein, protein2gene=protein2gene, configuration=config)

    # Building Translational Protein Sector
    translation_enzyme = TransEnzymeSector(id_list = ['R11'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=config)

    # Building Unused Enzyme Sector
    unused_enzyme = UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=config)

    # Building Membrane Sector
    alpha_numbers_dict = {"E1": 12,
                          "E2": 8,
                          "E3": 1,
                          "E4": 12,
                          "E5": 40,
                          "E7": 15,
                          "E8": 30}

    enzyme_location = {"E1": "Unknown",
                      "E2": "Cell membrane",
                      "E3": "Cytoplasm",
                      "E4": "Cytoplasm",
                      "E5": "Cell membrane",
                      "E7": "Cytoplasm",
                      "E8": "Cytoplasm"}

    membrane_sector = MembraneSector(area_avail_mu=0.1, area_avail_0=0.005, alpha_numbers_dict=alpha_numbers_dict,
                                     enzyme_location=enzyme_location, max_area=0.132815)

    # Building the toy_pam
    model = load_json_model('Models/toy_model_TS_v1.json')
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints',
                      sensitivity=sensitivity,
                      active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme,
                      membrane_sector=membrane_sector,
                      p_tot=0.2, configuration=config)

    return pamodel

### 1 Useful functions
def calculate_sensitivities(pamodel):
    glc_uptake_rates = np.linspace(1, 10, 20)
    Ccsc = []
    Cesc = []
    y_axis = []
    fluxes = []

    for glc in glc_uptake_rates:
        print('glucose uptake rate ', glc, ' mmol/gcdw/h')
        with pamodel:
            # change glucose uptake rate
            pamodel.change_reaction_bounds(rxn_id='R2',
                                           lower_bound=0, upper_bound=glc)
            # solve the model
            sol_pam = pamodel.optimize()
            fluxes.append(sol_pam.fluxes)
            if pamodel.solver.status == 'optimal': y_axis += [glc]
            # save data
            Ccsc_new = list()

            if pamodel.solver.status == 'optimal':
                capacity_coeff = pamodel.capacity_sensitivity_coefficients
                for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector', 'membrane']:
                    Ccsc_new += capacity_coeff[capacity_coeff['constraint'] == csc].coefficient.to_list()

                Ccsc += [Ccsc_new]

                enzyme_coeff = pamodel.enzyme_sensitivity_coefficients
                Cesc += [enzyme_coeff.coefficient.to_list()]

                print('Sum of capacity sensitivity coefficients: \t \t \t \t \t \t', round(sum(Ccsc_new), 6))
                print('Sum of enzyme sensitivity coefficients: \t \t \t \t \t \t', round(sum(Cesc[-1]), 6), '\n')

    return {'Ccsc': Ccsc, 'Cesc': Cesc, 'y_axis': y_axis, 'fluxes': fluxes, 'capacity coefficients': capacity_coeff,
            'enzyme coefficients': enzyme_coeff}


#
def make_heatmap_subfigure(results, csc_matrix, esc_matrix, x_csc, x_esc, yaxis, fig, grdspc,
                           ylabels=True, xlabels=False, cbar=True, title=None, fontsize=16,
                           vmin=-1.5, vmax=1.5, annotate=None, phenotype_data=None, cmap=None
                           # cmap = plt.cm.get_cmap('viridis')
                           ):
    # fig = plt.figure()
    if cmap is None:
        # Create separate colormaps for positive and negative values and a color for zero
        colors_neg = plt.cm.Blues(np.linspace(1, 0.3, 128))
        #         colors_pos = plt.cm.PuOr(np.linspace(0.5, 1, 128))#plt.cm.Reds(np.linspace(0, 0.5, 128))
        #         colors_pos = plt.cm.PuOr(np.linspace(0.5, 0, 64))  # Use part of the PuOr colormap
        #         colors_pos = np.vstack((colors_pos, plt.cm.PuOr(np.linspace(0.5, 1, 64))))
        #         colors_pos = lighten_colormap(plt.cm.plasma)(np.linspace(0, 0.5, 64))  # Use part of the PuOr colormap
        #         colors_pos = np.vstack((colors_pos, plt.cm.plasma(np.linspace(0.5, 1, 64))))
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

    esc_ax = fig.add_subplot(gs[2, 0])  # ESC heatmap
    acetate_ax = fig.add_subplot(gs[0, 0])  # acetate production
    csc_ax = fig.add_subplot(gs[1, 0], sharex=esc_ax)  # CSC heatmap
    if cbar:
        cbar_ax = fig.add_subplot(gs[1:, 1])  # colorbar

    # add annotation for subfigure (A or B)
    if annotate is not None:
        acetate_ax.annotate(annotate, xy=(2, 1), xycoords='data',
                            xytext=(-0.05, 1.5), textcoords='axes fraction',
                            va='top', ha='left', fontsize=fontsize * 1.5, weight='bold')

    glc_fluxes = [-sim.R2 for sim in results['fluxes']]

    # add arrow indicating growth regime
    # 0. remove the box to improve readability of the text
    acetate_ax.spines['top'].set_visible(False)
    acetate_ax.spines['right'].set_visible(False)

    # CSC heatmap
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

    return fig


#

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


#
# adjust labels for better readibility
def adjust_heatmap_labels(labels, pamodel):
    new_labels = labels.copy()

    for i, label in enumerate(labels):
        if 'EX_glc__D_e' in label or label[:-3] == 'EX_glc__D_e':
            if label[-1] == 'B':
                new_labels[i] = 'EX_glc_' + label[-2:]
            else:
                new_labels[i] = 'EX_glc_lb'
        if label == 'TotalProteinConstraint_proteome':
            new_labels[i] = 'Protein pool'
        if label[0].isdigit():  # all enzyme ids start with a digit
            rxn_ids = pamodel.get_reactions_with_enzyme_id(label)
            rxn_name = pamodel.reactions.get_by_id(rxn_ids[-1]).name.split('(')[0]
            new_labels[i] = '\n'.join([part for part in rxn_name.split(' ')])

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


#
def find_top5_sensitivities(Cv, x_axis, yaxis, threshold=0.005):
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

def parse_x_axis_heatmap(capacity_coeff, enzyme_coeff):
    x_axis_csc = []

    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector', 'membrane']:
        if csc == 'flux_ub' or csc == 'flux_lb':
            x_axis_csc += [coef + '_' + csc for coef in
                           capacity_coeff[capacity_coeff['constraint'] == csc].rxn_id.to_list()]
        else:
            x_axis_csc += [coef + '_' + csc for coef in capacity_coeff[
                capacity_coeff['constraint'] == csc].enzyme_id.to_list()]

    x_axis_esc = enzyme_coeff.enzyme_id.to_list()
    return x_axis_csc, x_axis_esc

BIOMASS_RXNID = Config.BIOMASS_REACTION
DATA_DIR = 'Data'  # os.path.join(os.path.split(os.getcwd())[0], 'Data')
PAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'proteinAllocationModel_iML1515_EnzymaticData_py_new_parsed.xls')
glc_uptake_rates = list(np.linspace(0.5, 10, 20))


if __name__ == "__main__":
    ### PAM simulations
    pamodel = build_toy_model()
    pamodel.objective = 'R11'
    results_pam = calculate_sensitivities(pamodel)
    x_axis_csc_pam, x_axis_esc_pam = parse_x_axis_heatmap(results_pam['capacity coefficients'],
                                                          results_pam['enzyme coefficients'])

    # get nonzero sensitivities
    csc_nonzero_pam, x_csc_nonzero_pam = find_nonzero_sensitivities(results_pam['Ccsc'], x_axis=x_axis_csc_pam)
    esc_nonzero_pam, x_esc_nonzero_pam = find_nonzero_sensitivities(results_pam['Cesc'], x_axis=x_axis_esc_pam)
    csc_nonzero_pam_t = np.transpose(np.array(csc_nonzero_pam))
    esc_nonzero_pam_t = np.transpose(np.array(esc_nonzero_pam))
    #
    cesc = results_pam['Cesc']

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

    ### 4.2 Load manually created plot with mapped sensitivities
    #
    image_path = os.path.join('Figures', 'Figure2_mapped_sensitivities.png')

    sensitivities_mapped = np.asarray(Image.open(image_path))
    #
    # create 2 plots: supplements and main text
    fontsize = 28
    width = 50
    height = 25
    # select colormap
    cmap = None  # plt.cm.get_cmap('magma')

    # gridspec inside gridspec
    fig = plt.figure(layout='constrained')

    gs0 = gridspec.GridSpec(1, 1, figure=fig)
    gs_pam = gs0[0]

    # adjust labels for better readibility
    x_csc_label_pam = adjust_heatmap_labels(x_csc_nonzero_pam)
    x_esc_label_pam = adjust_heatmap_labels(x_esc_nonzero_pam)

    fig_pam = make_heatmap_subfigure(results=results_pam, csc_matrix=csc_nonzero_pam_t, esc_matrix=esc_nonzero_pam_t,
                                     ylabels=True, xlabels=True, x_csc=x_csc_label_pam, x_esc=x_esc_label_pam,
                                     yaxis=glc_uptake_rates, fig=fig, grdspc=gs_pam,
                                     annotate='A', phenotype_data=pt_data, fontsize=fontsize, cmap=cmap)
    # fig_pam.subplots_adjust(left=0.3)
    # set common x axis title
    # fig_pam.xlabel('Glucose uptake rate [$mmol_{glc}/g_{CDW}/h$]', fontsize = fontsize*1.25)


    # plt.plasma()
    # fig.subplots_adjust(left=0.5)
    # fig.set_figwidth(width)
    # fig.set_figheight(height)
    # fig.align_labels()

    plt.show()

    # fig.savefig('Figures/test.png', dpi=200, bbox_inches='tight')

