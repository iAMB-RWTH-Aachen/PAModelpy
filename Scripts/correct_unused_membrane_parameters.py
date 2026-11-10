import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from src.PAModelpy.utils.pam_generation import set_up_pam
from src.PAModelpy.MembraneSector import UnusedMembraneSector
from Scripts.create_pamodel_from_diagnostics_file import create_pamodel_from_diagnostics_file

def scan_membrane_parameters(
    ue_slopes,
    unused_membrane_fractions,
    phenotype_data_path: str,
    phenotype_data_sheet: str,
    model
):
    results = []

    for slope in ue_slopes:
        for frac in unused_membrane_fractions:
            with model:
                print(f'Changing the unused enzyme slope and unused membrane fraction to: {slope} mg/gDW/h and {frac*100}%')
                model.sectors.get_by_id('MembraneSector').unused_membrane_sector.change_unused_membrane_parameters(
                    ups_slope=slope, 
                    unused_membrane_fraction=frac
                                            )

                df = calculate_simulation_error(
                    phenotype_data_path,
                    phenotype_data_sheet,
                    model
                )

                sse = df["total_error"].sum()

                results.append({
                    "ups_slope": slope,
                    "unused_membrane_fraction": frac,
                    "sse": sse
                })

                print('ups_slope:', slope)
                print('unused_membrane_fraction:', frac)
                print('sse:', sse)

    results_df = pd.DataFrame(results)

    heatmap_df = results_df.pivot(
        index="ups_slope",
        columns="unused_membrane_fraction",
        values="sse")

    return heatmap_df

def calculate_simulation_error(phenotype_data_path:str, phenotype_data_sheet:str, model):
    if phenotype_data_path is None:
        raise ValueError('Please provide the path to phenotype data')

    else:
        phenotype_df = pd.read_excel(phenotype_data_path, sheet_name=phenotype_data_sheet)
        phenotype_df = phenotype_df.reset_index(drop=True)
        phenotype_df["exp_id"] = phenotype_df.index
        
        # Rename the columns to avoid naming problem
        phenotype_df = phenotype_df.rename(columns={
        "BIOMASS_Ec_iML1515_core_75p37M": "growth_exp",
        "EX_ac_e": "acetate_exp",
        "EX_o2_e": "o2_exp",
        "EX_co2_e": "co2_exp"
        })

        # Drop strain and reference columns
        phenotype_df = phenotype_df.drop(['Reference', 'Strain'], axis=1)

    glc_uptake_rates = list(phenotype_df['EX_glc__D_e'])

    results = []

    for i, glc in enumerate(glc_uptake_rates):
        with model:
            # change glucose uptake rate
            model.change_reaction_bounds(rxn_id='EX_glc__D_e', upper_bound=glc, lower_bound=glc)
            # disable pyruvate formate lyase (inhibited by oxygen)
            model.change_reaction_bounds(rxn_id='PFL', upper_bound=0)
            # solve the model
            sol_pam = model.optimize()
            # save data
            results.append({
                "exp_id": i,
                "EX_glc_D_e_sim": glc,
                "growth_sim": sol_pam.fluxes[model.BIOMASS_REACTION], 
                "acetate_sim": sol_pam.fluxes["EX_ac_e"],
                "o2_sim": sol_pam.fluxes["EX_o2_e"],
                "co2_sim": sol_pam.fluxes["EX_co2_e"]
            })
    
    sim_df = pd.DataFrame(results)

    # Merge simulation data with experimental phenotype data
    merged_df = sim_df.merge(phenotype_df, on="exp_id")
    merged_df_with_error = compute_sse(df=merged_df)

    return merged_df_with_error

def compute_sse(df, biomass_weight: float = 0, acetate_weight: float = 10):
    eps = 1e-8  # avoid division by zero

    # Characteristic scales
    n_data = len(df['growth_sim'])

    # Normalized squared errors
    df["growth_err"] = (df["growth_sim"] - df["growth_exp"])**2 / n_data

    df["acetate_err"] = (df["acetate_sim"] - df["acetate_exp"])**2 / n_data

    # Weighted total error
    df["total_error"] = (
        biomass_weight * df["growth_err"] +
        acetate_weight * df["acetate_err"]
    )

    return df

def plot_and_save_heatmap(
                            heatmap_df,
                            output_path="heatmap.png",
                            title="Membrane parameter scan",
                        ):
    """
    Plot and save a heatmap.

    Parameters
    ----------
    heatmap_df : pandas.DataFrame
        DataFrame returned from pivot(), where
        - rows = UPS slope
        - columns = unused membrane fraction
        - values = SSE

    save_path : str
        Path where the figure will be saved.

    title : str
        Figure title.

    value_format : str
        Format string used for displaying the values in each cell.
    """

    heatmap = heatmap_df.astype(float)

    fig, ax = plt.subplots(figsize=(12, 8))

    im = ax.imshow(
        heatmap.values,
        origin="lower",
        aspect="auto",
        cmap="viridis"
    )

    # Colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Total error")

    # Axis ticks
    ax.set_xticks(np.arange(len(heatmap.columns)))
    ax.set_xticklabels(np.round(heatmap.columns, 3))

    ax.set_yticks(np.arange(len(heatmap.index)))
    ax.set_yticklabels(np.round(heatmap.index, 3))

    ax.set_xlabel("Unused membrane fraction")
    ax.set_ylabel("UPS slope")
    ax.set_title(title)

    # Write the value inside every cell
    for row in range(heatmap.shape[0]):
        for col in range(heatmap.shape[1]):

            value = heatmap.iloc[row, col]

            if np.isnan(value):
                text = "NA"
            else:
                text = format_value(value)

            ax.text(
                col,
                row,
                text,
                ha="center",
                va="center",
                color="white",
                fontsize=6
            )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)

def format_value(v):
    'Reformat long values into the shortest possible form'
    if np.isnan(v):
        return "NA"
    if abs(v) >= 1e4:
        return f"{v:.0e}"   # no decimals, shortest form
    return f"{v:.2g}"

def find_best_parameter_combination(heatmap):
    min_index = np.unravel_index(
        np.nanargmin(heatmap.values),
        heatmap.shape
    )

    best_slope = heatmap.index[min_index[0]]
    best_fraction = heatmap.columns[min_index[1]]
    best_sse = float(heatmap.iloc[min_index])

    return best_slope, best_fraction, best_sse
    
if __name__ == "__main__":
    # Define paths
    phenotype_data_path = "Data/Ecoli_phenotypes/Ecoli_phenotypes_py_rev.xlsx"
    pam_info_file = 'Results/PAM_parametrizer/Diagnostics_files/2026_06_10/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
    model_path = 'Models/iML1515.xml'

    # Set up mcpam
    mcpam = set_up_pam(pam_info_file=pam_info_file, 
                       model=model_path,
                       sensitivity=False, 
                       membrane_sector=True,
                       separate_memprot_from_tpc=True,
                       total_protein=0.241,
                       usable_area_fraction=0.6174,
                       unused_membrane_sector=UnusedMembraneSector()
                       )

    diagnostic_info_file = 'Results/PAM_parametrizer/Diagnostics_files/2026_06_10/pam_parametrizer_diagnostics_mciML1515_4.xlsx'
    mcpam = create_pamodel_from_diagnostics_file(file_path=diagnostic_info_file,
                                                 model=mcpam,
                                                 sheet_name='Best_Individuals')

    # Calculate error and visualize in heat map
    ue_slopes = np.arange(0.001, 0.015, 0.001)
    unused_membrane_fractions = np.arange(0.01, 0.18, 0.01)
    heatmap_df = scan_membrane_parameters(ue_slopes=[0.0075, 0.013],
                                          unused_membrane_fractions=[0.08, 0.16], 
                                          phenotype_data_path=phenotype_data_path,
                                          phenotype_data_sheet='Yields',
                                          model=mcpam)
    
    best_slope, best_fraction, best_sse = find_best_parameter_combination(heatmap=heatmap_df)
    print('Scan parameter finished. ')
    print('Best slope: ', best_slope, 'mmol/gDW/h')
    print('Best fraction: ', best_fraction)
    print('With sse: ', best_sse)

    plot_and_save_heatmap(heatmap_df, output_path="Results/PAM_parametrizer/Error_heat_maps/2026_06_10/error_heatmap_diagnostics_mciML1515_4")
