import pandas as pd
import os

from pam_generation import set_up_toy_pam, set_up_ecolicore_pam, set_up_ecoli_pam, parse_coefficients
from toy_ec_pam import print_heatmap
from numeric_error_estimation_schemes_vsc import (first_central_numeric_vsc_optimizations, first_forward_numeric_vsc_optimizations,
                                                  fcc_numeric_vsc_optimizations, first_forward_numeric_vsc_calculation,
                                                  first_central_numeric_vsc_calculation, fcc_numeric_vsc_calculation)

import time
import resource

GLC_UPTAKE = 9.82 #mmol/gcdw/hpamodel_inc.add_enzymes([eGFP_enzyme])
RESULT_DIR = os.path.join(os.path.split(os.getcwd())[0], 'Results')
RESULT_RESOURCE_FILE = os.path.join(RESULT_DIR, 'toy-core-full_model_performance_vsc.xls')
RESULT_VSC_FILE = os.path.join(RESULT_DIR, 'toy-core-full_model_VSC.xls')

"""
Script for estimation of time and memory usage for running sensitivity analyses on a toy, core and genome-scale model.

Comparing allocation coefficient method with numerical calculation.
"""
#define dataframe to store the results
resource_time_dataframe = pd.DataFrame(columns=['Model', 'Method', 'Sum', 'Time [s]', 'Resources [MByte]'])

#set-up models
toy_model_pam = set_up_toy_pam()
ecolicore_pam = set_up_ecolicore_pam()
ecoli_pam = set_up_ecoli_pam()


#set glucose uptake rate in the ecoli models to 9.81 for reproducible results
ecolicore_pam.change_reaction_bounds(rxn_id = 'EX_glc__D_e',
                                        lower_bound = -GLC_UPTAKE, upper_bound = -GLC_UPTAKE)
ecoli_pam.change_reaction_bounds(rxn_id = 'EX_glc__D_e',
                                        lower_bound = -GLC_UPTAKE, upper_bound = -GLC_UPTAKE)

coeff_dict = {}
#run simulations
for model in [toy_model_pam, ecolicore_pam]:#,
    model.optimize()

    x_axis_vsc = list()
    for vsc in ['rxn', 'enzyme', 'sector']:
        if vsc == 'rxn':
            x_axis_vsc += model.variable_sensitivity_coefficients[model.variable_sensitivity_coefficients['constraint'] == vsc].rxn_id.to_list()
        else:
            x_axis_vsc += model.variable_sensitivity_coefficients[
                model.variable_sensitivity_coefficients['constraint'] == vsc].enzyme_id.to_list()

    df_coeff = pd.DataFrame(columns=['method']+x_axis_vsc)
    #fcc
    time_start = time.perf_counter()
    fcc_vsc = fcc_numeric_vsc_optimizations(model)
    time_elapsed = (time.perf_counter() - time_start)
    res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
    accuracy = sum(fcc_vsc)
    resource_time_dataframe.loc[len(resource_time_dataframe)] = [model.id, 'fcc',
                                                                 accuracy, time_elapsed, res_used_mb]
    df_coeff.loc[len(df_coeff)] = ['FCC']+fcc_vsc
    print(f'Done with calculatibng Flux Control Coefficients for the {model}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}')


    # #forward
    # time_start = time.perf_counter()
    # ffn_vsc = first_forward_numeric_vsc_optimizations(model)
    # time_elapsed = (time.perf_counter() - time_start)
    # res_used_mb =resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0/1024.0
    # accuracy = 1/sum(ffn_vsc)
    # resource_time_dataframe.loc[len(resource_time_dataframe)] = [model.id, 'first order forward finite difference',
    #                                                              accuracy, time_elapsed, res_used_mb]
    # df_coeff.loc[len(df_coeff)] = ['Forward'] + fcc_vsc
    # print(f'Done with calculating first order forward finite difference for the {model}.\nTime:\t\t {time_elapsed} seconds \nAccuracy:\t\t{accuracy}')

    #central
    time_start = time.perf_counter()
    model.optimize()
    fcn_vsc = first_central_numeric_vsc_optimizations(model)
    time_elapsed = (time.perf_counter() - time_start)
    res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
    accuracy = sum(fcn_vsc)
    resource_time_dataframe.loc[len(resource_time_dataframe)] = [model.id, 'first order central finite difference',
                                                                 accuracy, time_elapsed, res_used_mb]
    df_coeff.loc[len(df_coeff)] = ['central'] + fcc_vsc
    print(f'Done with calculating first order central finite difference for the {model}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}')


    # allocation coefficients
    time_start = time.perf_counter()
    model.optimize()
    Ccsc,Cvsc = parse_coefficients(model)
    time_elapsed = (time.perf_counter() - time_start)
    res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
    accuracy = sum(Cvsc)
    resource_time_dataframe.loc[len(resource_time_dataframe)] = [model.id, 'allocation coefficients',
                                                                 accuracy,time_elapsed, res_used_mb]
    df_coeff.loc[len(df_coeff)] = ['allocation'] + fcc_vsc
    print(f'Done with calculating Variable Allocation Coefficients for {model}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}\n')

    coeff_dict[model.id] = df_coeff

    try:
        print_heatmap(x_axis_vsc, [fcc_vsc, fcn_vsc, Cvsc], yaxis=['fcc', 'central', 'allocation'])
    except:
        continue

#due to numeric complexity, the performance of the simulations with the full scale model will be estimated based on a single simulation
dummy_rxn = model.reactions[0]
nmbr_ineq_constraints = len(ecoli_pam.reactions)*2+len(ecoli_pam.enzymes)*2

#fcc
time_start = time.perf_counter()
ecoli_pam.optimize()
obj = ecoli_pam.objective.value
fcc_fac = fcc_numeric_vsc_calculation(ecoli_pam, dummy_rxn, f'{dummy_rxn.id}_ub', obj)
time_elapsed = (time.perf_counter() - time_start)*nmbr_ineq_constraints
res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
resource_time_dataframe.loc[len(resource_time_dataframe)] = [ecoli_pam.id, 'fcc',
                                                                 '', time_elapsed, res_used_mb]
print(f'Done with calculating Flux Control Coefficients for the {ecoli_pam}.\nEstimated time:\t\t {time_elapsed} seconds\n')

# #forward
# time_start = time.perf_counter()
# ecoli_pam.optimize()
# ffn_fac = first_forward_numeric_fac_calculation(ecoli_pam, dummy_rxn, f'{dummy_rxn.id}_ub', obj)
# time_elapsed = (time.perf_counter() - time_start)*nmbr_ineq_constraints
# res_used_mb =resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0/1024.0
# resource_time_dataframe.loc[len(resource_time_dataframe)] = [ecoli_pam.id, 'first order forward finite difference',
#                                                                  '', time_elapsed, res_used_mb]
# print(f'Done with calculating first order forward finite difference for the {ecoli_pam}.\nEstimated time:\t\t {time_elapsed} seconds\n')

#central
time_start = time.perf_counter()
ecoli_pam.optimize()
fcn_fac = first_central_numeric_vsc_calculation(ecoli_pam, dummy_rxn, f'{dummy_rxn.id}_ub', obj)
time_elapsed = (time.perf_counter() - time_start)*nmbr_ineq_constraints*2
res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
resource_time_dataframe.loc[len(resource_time_dataframe)] = [ecoli_pam.id, 'first order central finite difference',
                                                                 '', time_elapsed, res_used_mb]
print(f'Done with calculating first order central finite difference for the {ecoli_pam}.\nEstimated time:\t\t {time_elapsed} seconds\n')

# allocation coefficients
time_start = time.perf_counter()
ecoli_pam.optimize()
Ccac,Cfac = parse_coefficients(ecoli_pam)
time_elapsed = (time.perf_counter() - time_start)
res_used_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0
accuracy = sum(Cfac)
resource_time_dataframe.loc[len(resource_time_dataframe)] = [ecoli_pam.id, 'allocation coefficients',
                                                                 accuracy,time_elapsed, res_used_mb]
print(f'Done with calculating Variable Allocation Coefficients for {ecoli_pam}.\nTime:\t\t {time_elapsed} seconds \nSum:\t\t{accuracy}\n')

for row in resource_time_dataframe.iterrows():
    print(row)



#save the performance data
#check if there already is a Results directory
#otherwise make a results directory
if not os.path.exists(RESULT_DIR):
    os.mkdir(RESULT_DIR)
resource_time_dataframe.to_excel(RESULT_RESOURCE_FILE, engine= 'openpyxl')

writer = pd.ExcelWriter('MyFile.xlsx', engine='openpyxl')

#save each coeff matrix on a specific sheet
for sheet, frame in  coeff_dict.items():
    frame.to_excel(writer, sheet_name = sheet)