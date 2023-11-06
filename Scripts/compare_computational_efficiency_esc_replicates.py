import pandas as pd
import numpy as np
import os
import openpyxl

from pam_generation import set_up_toy_pam, set_up_ecolicore_pam, set_up_ecoli_pam, parse_esc
from numeric_error_estimation_schemes_esc import (first_central_numeric_esc_optimizations,
                                                  fcc_numeric_esc_optimizations,
                                                  first_central_numeric_esc_calculation, fcc_numeric_esc_calculation)

import time
import resource

NUM_REPLICATES = 5
GLC_UPTAKE = 9.82 #mmol/gcdw/hpamodel_inc.add_enzymes([eGFP_enzyme])
RESULT_DIR = os.path.join(os.path.split(os.getcwd())[0], 'Results')
RESULT_RESOURCE_FILE = os.path.join(RESULT_DIR, 'toy-core-full_model_performance_esc.xlsx')
RESULT_VSC_FILE = os.path.join(RESULT_DIR, 'toy-core-full_model_ESC.xlsx')

"""
Script for estimation of time and memory usage for running sensitivity analyses on a toy, core and genome-scale model.

Comparing allocation coefficient method with numerical calculation.
"""
#define dataframe to store the results
resource_time_dataframe = pd.DataFrame(columns=['Model', 'Method', 'Sum', 'Time [s]', 'sd'])

resource_time_dataframe_fcc = pd.DataFrame(columns=[ 'iteration','Model', 'Sum', 'Time [s]'])
resource_time_dataframe_fcn = pd.DataFrame(columns=[ 'iteration','Model', 'Sum', 'Time [s]'])
resource_time_dataframe_esc = pd.DataFrame(columns=['iteration', 'Model', 'Sum', 'Time [s]'])


#set-up models
toy_model_pam = set_up_toy_pam()
ecolicore_pam = set_up_ecolicore_pam()
ecoli_pam = set_up_ecoli_pam()


#set glucose uptake rate in the ecoli models to 9.81 for reproducible results
toy_model_pam.change_reaction_bounds(rxn_id = 'R1',
                                        lower_bound = 0, upper_bound = 0.1)
ecolicore_pam.change_reaction_bounds(rxn_id = 'EX_glc__D_e',
                                        lower_bound = -GLC_UPTAKE, upper_bound = -GLC_UPTAKE)
ecoli_pam.change_reaction_bounds(rxn_id = 'EX_glc__D_e',
                                        lower_bound = -GLC_UPTAKE, upper_bound = -GLC_UPTAKE)

coeff_dict = {}
#run simulations
for model in [toy_model_pam, ecolicore_pam]:
    model.optimize()

    x_axis_esc = model.enzyme_sensitivity_coefficients.enzyme_id.to_list()

    df_coeff_fcc = pd.DataFrame(columns=['iteration'] + x_axis_esc)
    df_coeff_fcn = pd.DataFrame(columns=['iteration'] + x_axis_esc)
    df_coeff_esc = pd.DataFrame(columns=['iteration'] + x_axis_esc)

    for i in range(NUM_REPLICATES):
        #fcc
        time_start = time.perf_counter()
        fcc_esc = fcc_numeric_esc_optimizations(model)
        time_elapsed = (time.perf_counter() - time_start)
        accuracy = sum(fcc_esc)
        resource_time_dataframe_fcc.loc[len(resource_time_dataframe_fcc)] = [i, model.id, accuracy, time_elapsed]
        df_coeff_fcc.loc[len(df_coeff_fcc)] = [i] + fcc_esc
        print(f'Done with calculatibng Flux Control Coefficients for the {model}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}')


        #central
        time_start = time.perf_counter()
        model.optimize()
        fcn_esc = first_central_numeric_esc_optimizations(model)
        time_elapsed = (time.perf_counter() - time_start)
        accuracy = sum(fcn_esc)
        resource_time_dataframe_fcn.loc[len(resource_time_dataframe_fcn)] = [i, model.id, accuracy, time_elapsed]
        df_coeff_fcn.loc[len(df_coeff_fcn)] = [i] + fcn_esc
        print(f'Done with calculating first order central finite difference for the {model}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}')


        # allocation coefficients
        time_start = time.perf_counter()
        model.optimize()
        Cesc = parse_esc(model)
        time_elapsed = (time.perf_counter() - time_start)
        res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
        accuracy = sum(Cesc)
        resource_time_dataframe_esc.loc[len(resource_time_dataframe_esc)] = [i, model.id, accuracy, time_elapsed]
        df_coeff_esc.loc[len(df_coeff_esc)] = [i] + Cesc
        print(f'Done with calculating Enzyme Sensitivity Coefficients for {model}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}\n')

    coeff_dict[model.id] = {'fcc':df_coeff_fcc, 'fcn':df_coeff_fcn, 'evsc':df_coeff_esc}


#due to numeric complexity, the performance of the simulations with the full scale model will be estimated based on a single simulation
dummy_enz = ecoli_pam.enzyme_variables[0]
nmbr_ineq_enz_constraints = len(ecoli_pam.enzymes) * 2

for i in range(NUM_REPLICATES):
    #fcc
    time_start = time.perf_counter()
    ecoli_pam.optimize()
    obj = ecoli_pam.objective.value
    fcc_esc = fcc_numeric_esc_calculation(ecoli_pam, dummy_enz, f'{dummy_enz.id}_max', obj)
    time_elapsed = (time.perf_counter() - time_start) * nmbr_ineq_enz_constraints
    resource_time_dataframe_fcc.loc[len(resource_time_dataframe_fcc)] = [i, ecoli_pam.id, '',time_elapsed]
    print(f'Done with calculating Flux Control Coefficients for the {ecoli_pam}.\nEstimated time:\t\t {time_elapsed} seconds\n')

    #central
    time_start = time.perf_counter()
    ecoli_pam.optimize()
    fcn_esc = first_central_numeric_esc_calculation(ecoli_pam, dummy_enz, f'{dummy_enz.id}_max', obj)
    time_elapsed = (time.perf_counter() - time_start) * nmbr_ineq_enz_constraints * 2
    resource_time_dataframe_fcn.loc[len(resource_time_dataframe_fcn)] = [i, ecoli_pam.id, '', time_elapsed]
    print(f'Done with calculating first order central finite difference for the {ecoli_pam}.\nEstimated time:\t\t {time_elapsed} seconds\n')

    # allocation coefficients
    time_start = time.perf_counter()
    ecoli_pam.optimize()
    Cesc = parse_esc(ecoli_pam)
    time_elapsed = (time.perf_counter() - time_start)
    accuracy = sum(Cesc)
    resource_time_dataframe_esc.loc[len(resource_time_dataframe_esc)] = [i, ecoli_pam.id, accuracy, time_elapsed]
    print(f'Done with calculating Enzyme Sensitivity Coefficients for {ecoli_pam}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}\n')

#calculate averages and standard deviation
# resource_time_dataframe = pd.DataFrame(columns=['Model', 'Method', 'Sum', 'Time [s]', 'sd'])
for resource in zip([resource_time_dataframe_fcc, resource_time_dataframe_fcn, resource_time_dataframe_esc], ['fcc', 'fcn', 'sEnz']):
    for model_id in [toy_model_pam.id, ecolicore_pam.id, ecoli_pam.id]:
        resource_df = resource[0]
        resource_info = resource_df[resource_df['Model'] == model_id]
        method = resource[1]
        #cannot get a mean of the sum for the full model for fcc and fcn, because sums are not calculated
        if method != 'sEnz' and model_id == ecoli_pam.id:
            sum = ''
        else:
            sum = resource_info['Sum'].mean()
        time_avg = resource_info['Time [s]'].mean()
        time_sd = resource_info['Time [s]'].std()
        resource_time_dataframe.loc[len(resource_time_dataframe)] = [model_id, method, sum, time_avg, time_sd]


#save the performance data
#check if there already is a Results directory
#otherwise make a results directory
if not os.path.exists(RESULT_DIR):
    os.mkdir(RESULT_DIR)

#saving the resource information
resource_time_dataframe.to_excel(RESULT_RESOURCE_FILE, sheet_name='results', index=False, engine='openpyxl')

with pd.ExcelWriter(RESULT_RESOURCE_FILE, engine='openpyxl', mode='a') as writer:
    # Load the existing Excel file
    writer.book = openpyxl.load_workbook(RESULT_RESOURCE_FILE)

    # Add DataFrames for Method 1 and Method 2 to the Excel file
    resource_time_dataframe_esc.to_excel(writer, sheet_name='sEnz', index=False)
    resource_time_dataframe_fcc.to_excel(writer, sheet_name='flux_control_coefficient', index=False)
    resource_time_dataframe_fcn.to_excel(writer, sheet_name='first_central_difference', index=False)

#saving the enzyme VSC information
coeff_dict[ecolicore_pam.id]['fcc'].to_excel(RESULT_VSC_FILE, sheet_name='flux_control_coefficient', index=False)

with pd.ExcelWriter(RESULT_VSC_FILE, engine='openpyxl', mode='a') as writer:
    # Load the existing Excel file
    writer.book = openpyxl.load_workbook(RESULT_VSC_FILE)

    # Add DataFrames for Method 1 and Method 2 to the Excel file
    coeff_dict[ecolicore_pam.id]['fcn'].to_excel(writer, sheet_name='first_central_difference', index=False)
    coeff_dict[ecolicore_pam.id]['evsc'].to_excel(writer, sheet_name='sEnz', index=False)