from cobra import Model
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
from Scripts.create_pamodel_from_diagnostics_file import change_prot_kcats

def change_set_of_kcats_using_excel_sheet(models:list, prot_file_path:str, sheet:str):
     prot_df = pd.read_excel(io = prot_file_path, sheet_name = sheet)
     for model in models:
        model = change_prot_kcats(prot_df=prot_df, model=model)

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
                    # model.reactions.EX_glc__D_e.lower_bound = -glc
                    model.change_reaction_bounds(rxn_id = 'EX_glc__D_e', 
                                            lower_bound = -glc, upper_bound = -glc)
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
                    # model.reactions.EX_glc__D_e.lower_bound = -glc
                    model.change_reaction_bounds(rxn_id = 'EX_glc__D_e', 
                                            lower_bound = -glc, upper_bound = -glc)
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

    return fig

def run_simulations_pam_mcpam_w_different_areas(models, print_area:bool=False, type:str="full scale"):
    fontsize = 25
    labelsize = 15

    # load phenotype data from excel file
    pt_data = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes','Ecoli_phenotypes_py_rev.xls'),
                            sheet_name='Yields', index_col=None)

    # Define the biomass name based on the used model
    if type == "full scale":
        biomass_name = 'BIOMASS_Ec_iML1515_core_75p37M'
        max_area_list = np.linspace(0.01, 0.1, 10)
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
                    # model.reactions.EX_glc__D_e.lower_bound = -glc
                    model.change_reaction_bounds(rxn_id = 'EX_glc__D_e', 
                                            lower_bound = -glc, upper_bound = -glc)
                    # disable pyruvate formate lyase (inhibited by oxygen)
                    model.reactions.PFL.upper_bound = 0
                    # solve the model
                    # model.objective = 'EX_ac_e'
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

                print(model.constraints.membrane)
                occupied_area, available_area = model.sectors.get_by_id(
                    'MembraneSector').calculate_occupied_membrane(model)
                print('Available area: ', available_area, 'um2')
                print('Occupied area: ', occupied_area, 'um2')
                print('Occupied area: ', occupied_area / available_area * 100, '%')
                print('Growth rate: ', model.objective.value)
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

def get_info_for_proteins(mcpam, pam_info_path, protein_info_path) -> None:
    '''
    Create an excel sheet with following protein (enzyme) information:
        - enzyme_id
        - reaction_id
        - forward kcat
        - backward kcat
        - Alpha number
        - Occupied area
        - Contribution to protein pool

    Return:
        None, create an excel sheet in the provided data path
    '''
    file_name = os.path.basename(pam_info_path)  # 'proteinAllocationModel_EnzymaticData_iML1515_10.xlsx'
    file_stem = os.path.splitext(file_name)[0]  # 'proteinAllocationModel_EnzymaticData_iML1515_10'
    match = re.search(r'_(\d+)$', file_stem)
    number = match.group(1) if match else None
    prot_occupancy_df = mcpam.sectors.get_by_id('MembraneSector').calculate_occupied_membrane(mcpam, get_df=True)

    # Write excel datasheet
    with pd.ExcelWriter(protein_info_path, engine='openpyxl', mode='a') as writer:
        # Write the new DataFrame to a new sheet
        prot_occupancy_df.to_excel(writer, sheet_name=f'enzymatic_file_{number}', index=True)

def get_missing_backward_kcats(pamodel):

    missing_kcats_b = []

    for rxn in pamodel.reactions:
        reversibility = rxn.reversibility

        if reversibility:
            enzymes = pamodel.get_enzymes_with_reaction_id(rxn.id)

            if enzymes is not None:
                for enzyme in enzymes:
                    for rxn_id, kcats in enzyme.enzyme_variable.kcats.items():
                        if kcats['b'] == 0:
                            missing_kcats_b.append({
                                'enzyme_id': enzyme.id,
                                'Reaction': rxn_id,
                                'kcat_f': kcats.get('f'),
                                'kcat_b': kcats.get('b')
                            })
    df = pd.DataFrame(missing_kcats_b)
    
    return df

def fill_missing_backward_kcats(df:pd.DataFrame)->pd.DataFrame:
    for i, row in df.iterrows():
        df['kcat_b'].iloc[i] = df['kcat_f'].iloc[i]

    return df
    
