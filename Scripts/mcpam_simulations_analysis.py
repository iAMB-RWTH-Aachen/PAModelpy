import cobra
import pandas as pd
import numpy as np
import os
from typing import Union
import seaborn as sns
import json
from matplotlib import pyplot as plt
from cobra.flux_analysis import single_gene_deletion

# load mcPAMpy modules
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.PAModelpy.MembraneSector import MembraneSector
from src.PAModelpy.configuration import Config
from Scripts.mcpam_generation_uniprot_id import (parse_reaction2protein,
                                                 set_up_ecolicore_pam, set_up_ecolicore_mcpam,
                                                 set_up_ecolicore_mcpam_new_surface_parameter,
                                                 set_up_ecoli_pam, set_up_ecoli_mcpam)

def compare_mu_for_different_sensitivities_ecolicore_pam():
    rxn2protein_new = {}
    config = Config()
    config.reset()
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]

    enzyme_info_path = "Data/mcPAM_iML1515_EnzymaticData.xlsx"
    active_enzyme_info = pd.read_excel(enzyme_info_path, sheet_name='mcPAM_data')

    model = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))
    rxn2protein, protein2gene = parse_reaction2protein(active_enzyme_info, model)

   # building the different sectors
    # building translational protein sector
    # translational protein sector parameter (substrate dependent)
    id_list_tps = ['EX_glc__D_e']
    tps_0 = [0.04992]  # g/gDW
    tps_mu = [-0.002944]  # g h/gDW -> transformed to match glucose uptake variable
    molmass_tps = [405903.94]  # g/mol

    # translational protein sector
    translation_enzyme_sector = TransEnzymeSector(
        id_list=id_list_tps,
        tps_0=tps_0,
        tps_mu=tps_mu,
        mol_mass=molmass_tps,
    )

    # building unused protein sector
    config.BIOMASS_REACTION = BIOMASS_REACTION
    id_list_ups = [BIOMASS_REACTION]
    ups_0 = [0.0407]  # g/gDW
    ups_mu = [-0.0214]  # g h/gDW -> negative relation with growth rate
    molmass_ups = [405903.94]  # g/mol

    unused_enzyme_sector = UnusedEnzymeSector(
        id_list=id_list_ups,
        ups_0=ups_0,
        ups_mu=ups_mu,
        mol_mass=molmass_ups,
    )

    # building membrane sector #addition
    membrane_info = pd.read_excel(enzyme_info_path, sheet_name='Membrane')

    area_avail_0 = membrane_info[membrane_info.Parameter == 'area_avail_0'].loc[1, 'Value']
    area_avail_mu = membrane_info[membrane_info.Parameter == 'area_avail_mu'].loc[2, 'Value']
    alpha_numbers_dict = active_enzyme_info.set_index(keys='uniprotID').loc[:, 'alpha_numbers'].to_dict()
    enzyme_location = active_enzyme_info.set_index(keys='uniprotID').loc[:, 'Location'].to_dict()

    membrane_sector = MembraneSector(area_avail_0=area_avail_0,
                                     area_avail_mu=area_avail_mu,
                                     alpha_numbers_dict=alpha_numbers_dict,
                                     enzyme_location=enzyme_location)

    # creating models by adding one reaction at a time
    df_mu = pd.DataFrame(columns=['True', 'False'])
    mu_true = []
    mu_false = []
    loop_count = 0

    for reaction, proteins in rxn2protein.items():
        model1 = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))
        model2 = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))

        for protein, value in proteins.items():
            if reaction not in rxn2protein_new:
                rxn2protein_new[reaction] = {}  # Initialize the nested dictionary for each reaction

                rxn2protein_new[reaction][protein] = value  # Add the protein-value pair to the nested dictionary
            else:
                rxn2protein_new[reaction][protein] = value
        active_enzyme_sector_1 = ActiveEnzymeSector(rxn2protein=rxn2protein_new, protein2gene=protein2gene,
                                                    configuration=config)

        active_enzyme_sector_2 = ActiveEnzymeSector(rxn2protein=rxn2protein_new.copy(),
                                                    protein2gene=protein2gene.copy(),
                                                    configuration=config)

        pamodel1 = PAModel(id_or_model=model1, p_tot=TOTAL_PROTEIN_CONCENTRATION,
                           active_sector=active_enzyme_sector_1, translational_sector=translation_enzyme_sector,
                           unused_sector=unused_enzyme_sector, sensitivity=True, configuration=config,
                           membrane_sector=membrane_sector
                           )
        pamodel2 = PAModel(id_or_model=model2, p_tot=TOTAL_PROTEIN_CONCENTRATION,
                           active_sector=active_enzyme_sector_2, translational_sector=translation_enzyme_sector,
                           unused_sector=unused_enzyme_sector, sensitivity=False, configuration=config,
                           membrane_sector=membrane_sector
                           )

        pamodel1.objective = pamodel1.BIOMASS_REACTION
        pamodel2.objective = pamodel2.BIOMASS_REACTION

        pamodel1.optimize()
        pamodel2.optimize()

        obj_true = pamodel1.objective.value
        obj_false = pamodel2.objective.value

        if pamodel1.solver.status == 'optimal':
            obj_true = round(obj_true, 2)
            obj_false = round(obj_false, 2)

        mu_true.append(obj_true)
        mu_false.append(obj_false)

        if obj_true == obj_false:
            print(obj_true, '=', obj_false)

        # elif loop_count == 7:
        #     break
        else:
            print('Objective values are not the same')
            print('Stop at reaction', reaction, rxn2protein_new[reaction])

        loop_count += 1

    df_mu['True'] = mu_true
    df_mu['False'] = mu_false

    return df_mu

def change_kcats_for_multiple_enzyme_sets(models:list, enzyme_sets:list):
    for model in models:
        for enzyme_set in enzyme_sets.values():
            change_kcats_for_an_enzyme_set(model, enzyme_set)

def change_kcats_for_an_enzyme_set(model:PAModel, enzyme_set:dict):
    for enzyme, kcats in enzyme_set.items():
        model.change_kcat_value(enzyme, kcats)

def build_ecolicore_pam_and_mcpam(sensitivity:bool=True, max_area: float = 0.1):
    pam_core = set_up_ecolicore_pam(sensitivity=sensitivity)
    mcpam_core = set_up_ecolicore_mcpam(sensitivity=sensitivity, max_area=max_area)

    return pam_core, mcpam_core

def build_ecoli_pam_and_mcpam(sensitivity:bool=True, max_area: float = 0.1):
    pam_core = set_up_ecoli_pam(sensitivity=sensitivity)
    mcpam_core = set_up_ecoli_mcpam(sensitivity=sensitivity)

    return pam_core, mcpam_core

def set_objective_pam_mcpam(models:list, objective:str):
    for model in models:
        model.objective = objective

def optimize_pam_mcpam(models: list):
    for model in models:
        model.optimize()

def run_pam_mcpam_core_with_optimized_kcats(sensitivity:bool=True,
                                            enzyme_sets_name:str='enzyme_sets.json',
                                            print_area:bool=False,
                                            type:str='full scale'):

    # Resetting the configurations
    config = Config()
    config.reset()

    # Building the models
    if type == "full scale":
        pam, mcpam =build_ecoli_pam_and_mcpam(sensitivity=sensitivity)
    else:
        pam, mcpam = build_ecolicore_pam_and_mcpam(sensitivity=sensitivity)

    # List of kcat values to be changed
    enzyme_path = os.path.join('Data', enzyme_sets_name)
    with open(enzyme_path, 'r') as json_file:
        enzyme_sets = json.load(json_file)

    # Changing the kcats
    models = [pam, mcpam]
    change_kcats_for_multiple_enzyme_sets(models, enzyme_sets)

    # Setting the objectives and Optimizing the models
    if type == "full scale":
        BIOMASS_REACTION = 'BIOMASS_Ec_iML1515_core_75p37M'
    else:
        BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'

    set_objective_pam_mcpam(models, BIOMASS_REACTION)
    optimize_pam_mcpam(models)

    ## Printing the objective values
    for model, config in zip(models, ['pam', 'mcpam']):
        print(f'{config} objective value:', model.objective.value)

    #Printing the occupied area for mcpam core
    if print_area:
        occupied_area, available_area = mcpam.calculate_occupied_membrane()
        print(f'occupied area {occupied_area / available_area * 100}%')

    return pam, mcpam

def run_simulation_pam_mcpam(models, type:str="full scale"):
    fontsize = 25
    labelsize = 15

    # load phenotype data from excel file
    pt_data = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls'),
                            sheet_name='Yields', index_col=None)

    # Define the biomass name based on the used model
    if type == "full scale":
        biomass_name = 'BIOMASS_Ec_iML1515_core_75p37M'
    else:
        biomass_name = 'BIOMASS_Ecoli_core_w_GAM'

    # extract reaction specific data
    rxn_to_pt = {}
    rxn_transform = {
        'EX_ac_e': 'EX_ac_e',
        'EX_co2_e': 'EX_co2_e',
        'EX_o2_e': 'EX_o2_e',
        biomass_name: 'BIOMASS_Ec_iML1515_core_75p37M'
    }
    for rxn_id, pt_id in rxn_transform.items():
        rxn_to_pt[rxn_id] = pt_data[['EX_glc__D_e', pt_id]].dropna().rename(columns={pt_id: rxn_id})

    glc_uptake_rates = np.linspace(0.5, 14, 25)

    # Initializing fluxes and concentrations for pam and mcpam
    fluxes_dict = {}
    concentrations_dict = {}

    for model, config in zip(models, ["PAM", "mcPAM"]):

        # disable pyruvate formate lyase (inhibited by oxygen)
        model.change_reaction_bounds(rxn_id='PFL', upper_bound=0)

        fluxes_list = []
        concentrations_list = [0]

        if config == "PAM":  # simulating pam
            for glc in glc_uptake_rates:
                with model:
                    # change glucose uptake rate
                    model.reactions.EX_glc__D_e.lower_bound = -glc
                    # disable pyruvate formate lyase (inhibited by oxygen)
                    model.reactions.PFL.upper_bound = 0
                    # solve the model
                    sol_pam = model.optimize()
                    # save data
                    fluxes_list.append(sol_pam.fluxes)  # flux distributions
                    concentration = 0
                    for enz_var in model.enzyme_variables:
                        concentration += enz_var.concentration
                    concentrations_list.append(concentration)

            fluxes_dict[config] = fluxes_list
            concentrations_dict[config] = concentrations_list

        else:  # simulating mcpam
            for glc in glc_uptake_rates:
                with model:
                    # change glucose uptake rate
                    model.reactions.EX_glc__D_e.lower_bound = -glc
                    # disable pyruvate formate lyase (inhibited by oxygen)
                    model.reactions.PFL.upper_bound = 0
                    # solve the model
                    sol_pam = model.optimize()
                    # save data
                    fluxes_list.append(sol_pam.fluxes)  # flux distributions
                    concentration = 0
                    for enz_var in model.enzyme_variables:
                        concentration += enz_var.concentration
                    concentrations_list.append(concentration)

            fluxes_dict[config] = fluxes_list
            concentrations_dict[config] = concentrations_list


    # dictionary of colorblind friendly color palette
    sns.set_palette(("colorblind"))

    # plot flux changes with glucose uptake
    rxn_id = ['EX_ac_e', 'EX_co2_e', 'EX_o2_e', biomass_name]
    ax_title = {'EX_ac_e': 'Acetate Secretion',
                'EX_co2_e': 'CO2 Secretion',
                'EX_o2_e': 'O2 uptake',
                biomass_name: 'Biomass Production'}
    # rxn_id = ['EX_ac_e', 'EX_co2_e', 'EX_o2_e', 'BIOMASS_Ec_iML1515_core_75p37M']
    fig, axs = plt.subplots(2, 2, dpi=90)
    for r, ax in zip(rxn_id, axs.flatten()):
        # plot data
        if r in rxn_to_pt.keys():
            ax.scatter(abs(rxn_to_pt[r]['EX_glc__D_e']), abs(rxn_to_pt[r][r]),
                       color='firebrick', marker='o', s=30, linewidths=1.3,
                       facecolors=None, zorder=0,
                       label='Data')

        # plot simulation for pam core
        for model in fluxes_dict.keys():
            if model == 'PAM':
                ax.plot(glc_uptake_rates, [abs(f[r]) for f in fluxes_dict[model]],
                        label=f'{model}', linewidth=2.5, linestyle='--',
                        zorder=6)
            else:
                ax.plot(glc_uptake_rates, [abs(f[r]) for f in fluxes_dict[model]],
                        label=f'{model}', linewidth=2.5,
                        zorder=5)

        # options
        ax.tick_params(axis='both', which='major', labelsize=labelsize)
        ax.set_xlabel('glc uptake rate [mmol/gDW/h]', fontsize=fontsize)
        ax.set_ylabel('flux [mmol/gDW/h]', fontsize=fontsize * 0.8)
        ax.set_title(ax_title[r], fontsize=fontsize * 0.8, fontweight="bold")
        # set grid
        ax.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
        ax.set_axisbelow(True)
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.95), ncol=6, fontsize=fontsize * 0.65)

    # show legend
    fig.subplots_adjust(top=0.85, wspace=0.5, hspace=0.5)

    plt.show()

def run_simulations_pam_mcpam_w_different_areas(models, print_area:bool=False, type:str="full scale"):
    fontsize = 25
    labelsize = 15

    # load phenotype data from excel file
    pt_data = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes','Ecoli_phenotypes_py_rev.xls'),
                            sheet_name='Yields', index_col=None)

    # Define the biomass name based on the used model
    if type == "full scale":
        biomass_name = 'BIOMASS_Ec_iML1515_core_75p37M'
        max_area_list = np.linspace(0.1, 0.4, 4)
    else:
        biomass_name = 'BIOMASS_Ecoli_core_w_GAM'
        max_area_list = np.linspace(0.01, 0.04, 4)

    # extract reaction specific data
    rxn_to_pt = {}
    rxn_transform = {
        'EX_ac_e': 'EX_ac_e',
        'EX_co2_e': 'EX_co2_e',
        'EX_o2_e': 'EX_o2_e',
        biomass_name:'BIOMASS_Ec_iML1515_core_75p37M'
    }
    for rxn_id, pt_id in rxn_transform.items():
        rxn_to_pt[rxn_id] = pt_data[['EX_glc__D_e', pt_id]].dropna().rename(columns={pt_id: rxn_id})

    glc_uptake_rates = np.linspace(0.5, 14, 25)

    # Initializing fluxes and concentrations for pam and mcpam
    fluxes_dict = {}
    concentrations_dict = {}

    for model, config in zip(models, ["PAM", "mcPAM"]):

        # disable pyruvate formate lyase (inhibited by oxygen)
        model.change_reaction_bounds(rxn_id='PFL', upper_bound=0)

        fluxes_list = []
        concentrations_list = [0]

        if config == "PAM": # simulating pam_core
            for glc in glc_uptake_rates:
                with model:
                    # change glucose uptake rate
                    model.reactions.EX_glc__D_e.lower_bound = -glc
                    # disable pyruvate formate lyase (inhibited by oxygen)
                    model.reactions.PFL.upper_bound = 0
                    # solve the model
                    sol_pam = model.optimize()
                    # save data
                    fluxes_list.append(sol_pam.fluxes)  # flux distributions
                    concentration = 0
                    for enz_var in model.enzyme_variables:
                        concentration += enz_var.concentration
                    concentrations_list.append(concentration)

            fluxes_dict[config] = fluxes_list
            concentrations_dict[config] = concentrations_list

        else: #simulating mcpam_core
            for area in max_area_list:
                for glc in glc_uptake_rates:
                    with model:
                        # change glucose uptake rate
                        model.reactions.EX_glc__D_e.lower_bound = -glc
                        # disable pyruvate formate lyase (inhibited by oxygen)
                        model.reactions.PFL.upper_bound = 0
                        # change max available membrane area
                        model.sectors.get_by_id('MembraneSector').change_available_membrane_area(area, model)
                        # solve the model
                        sol_pam = model.optimize()
                        # save data
                        fluxes_list.append(sol_pam.fluxes)  # flux distributions
                        concentration = 0
                        for enz_var in model.enzyme_variables:
                            concentration += enz_var.concentration
                        concentrations_list.append(concentration)

                area = float("{:.2f}".format(area))*100
                key = f'{config} {str(area)} % area'
                fluxes_dict[key] = fluxes_list
                concentrations_dict[config+str(area)] = concentrations_list
                fluxes_list = []
                concentrations_list = [0]

    # dictionary of colorblind friendly color palette
    sns.set_palette(("colorblind"))

    # plot flux changes with glucose uptake
    rxn_id = ['EX_ac_e', 'EX_co2_e', 'EX_o2_e', biomass_name]
    ax_title = {'EX_ac_e': 'Acetate Secretion',
                'EX_co2_e': 'CO2 Secretion',
                'EX_o2_e': 'O2 uptake',
                biomass_name: 'Biomass Production'}
    # rxn_id = ['EX_ac_e', 'EX_co2_e', 'EX_o2_e', 'BIOMASS_Ec_iML1515_core_75p37M']
    fig, axs = plt.subplots(2, 2, dpi=90)
    for r, ax in zip(rxn_id, axs.flatten()):
        # plot data
        if r in rxn_to_pt.keys():
            ax.scatter(abs(rxn_to_pt[r]['EX_glc__D_e']), abs(rxn_to_pt[r][r]),
                       color='black', marker='o', s=30, linewidths=1.3,
                       facecolors=None, zorder=0,
                       label='Data')

        # plot simulation for pam core
        for model in fluxes_dict.keys():
            if model == 'PAM':
                ax.plot(glc_uptake_rates, [abs(f[r]) for f in fluxes_dict[model]],
                        label=f'{model}', linewidth=2.5, linestyle='--',
                        zorder=6)
            else:
                ax.plot(glc_uptake_rates, [abs(f[r]) for f in fluxes_dict[model]],
                        label=f'{model}', linewidth=2.5,
                        zorder=5)

        # options
        ax.tick_params(axis='both', which='major', labelsize=labelsize)
        ax.set_xlabel('glc uptake rate [mmol/gDW/h]', fontsize=fontsize)
        ax.set_ylabel('flux [mmol/gDW/h]', fontsize=fontsize*0.8)
        ax.set_title(ax_title[r], fontsize=fontsize*0.8, fontweight="bold")
        # set grid
        ax.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
        ax.set_axisbelow(True)
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.95), ncol=6, fontsize=fontsize*0.65)

    # show legend
    fig.subplots_adjust(top=0.85, wspace=0.5, hspace=0.5)

    plt.show()

def perform_single_gene_ko_for_all_genes(model):

    # Performing the gene knockouts
    model.optimize()
    mu_ori = model.objective.value
    print("standard growth: ", mu_ori)

    deletion_results = single_gene_deletion(model)
    deletion_results['ids'] = deletion_results['ids'].astype(str).str.replace(r"[{''}]", '', regex=True)

    # Add gene names to the data frame
    data2_path = os.path.join('Data/proteinAllocationModel_iML1515_EnzymaticData_py_uniprot.xlsx')
    data2 = pd.read_excel(data2_path, sheet_name="ActiveEnzymes")

    data2 = data2[['m_gene', 'gene_name']].drop_duplicates()

    merged_data = pd.merge(deletion_results, data2, how="left", left_on="ids", right_on="m_gene")
    knockout_results = merged_data[['m_gene', 'gene_name', 'growth', 'status']]

    return knockout_results


def perform_and_plot_single_KO(model, genes_to_be_ko:list):

    # Determining the fluxes to be studied
    fluxes_name = ['PGI', 'PGL', 'FBA', 'TPI', 'GAPD', 'PGK', 'PYK', 'CS', 'ICDHyr', 'FUM', 'MDH', 'ICL', 'ACKr',
                   'LDH_D', 'G6PDH2r', 'EX_ac_e']
    fluxes = {}

    for gene in genes_to_be_ko:
        fluxes[gene] = {}

        for flux in fluxes_name:
            fluxes[gene][flux] = {}

    # Determining the fluxes for wild type
    fluxes['wild type'] = {}
    sol_pam = model.optimize()
    print("Wild type growth: ", model.objective.value)
    for flux in fluxes_name:
        fluxes['wild type'][flux] = sol_pam.fluxes[flux]

    # Performing KOs for every gene in the list
    for gene in genes_to_be_ko:
        ko_gene = getattr(model.genes, gene)
        ko_gene.knock_out()
        sol_pam = model.optimize()
        print(f'Strain {gene} growth: ', model.objective.value)

        for flux in fluxes_name:
            fluxes[gene][flux] = sol_pam.fluxes[flux]

    # Plotting the fluxes
    num_strains = len(fluxes)
    bar_width = 0.8 / num_strains  # Adjust width to fit all bars side by side
    x = np.arange(len(next(iter(fluxes.values()))))  # x positions for each flux type

    # Plot each strain's data
    for i, (strain, flux) in enumerate(fluxes.items()):
        plt.bar(
            x + i * bar_width,  # Shift each strain's bars by i * bar_width
            list(flux.values()),  # Heights of bars
            width=bar_width,  # Width of each bar
            align='center',
            alpha=0.5,
            label=strain
        )
    # Add labels and legend
    plt.xticks(x + bar_width * (num_strains - 1) / 2, list(next(iter(fluxes.values())).keys()))
    plt.xlabel('Fluxes')
    plt.ylabel('Values')
    plt.title('Flux Comparison')
    plt.legend()

    plt.show()

def get_memprot_data_in_mcpam(memprot_dict: dict):
    df_list = []
    for protein_group, (flux_dict, alpha_number) in memprot_dict.items():
        for reaction, flux_values in flux_dict.items():
            df_list.append({
                'Protein Group': protein_group,
                'Reaction': reaction,
                'Forward Flux': flux_values['f'],
                'Backward Flux': flux_values['b'],
                'Alpha Number': alpha_number
            })

    df = pd.DataFrame(df_list)

    return df