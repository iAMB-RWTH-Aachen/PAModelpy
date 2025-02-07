import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data1_path = os.path.join('Data/mcPAM_iML1515_EnzymaticData.xlsx')
data1 = pd.read_excel(data1_path, sheet_name='mcPAM_data_core')
core_data_path = os.path.join('Data/proteinAllocationModel_mc-core_EnzymaticData_241209_multi.xlsx')

# Unpivot the DataFrame
columns = ['uniprotID', 'rxnID', 'molMass', 'm_gene_reaction_rule', 'm_gene', 'direction', 'kcat']
data1 = data1.melt(
    id_vars=columns[0:5],                # Columns to keep as is
    value_vars=['kcat_f', 'kcat_b'],  # Columns to unpivot
    var_name='direction',             # Name of the new column for direction
    value_name='kcat'                 # Name of the new column for kcat values
)
data1["direction"] = data1["direction"].replace({"kcat_f": "f", "kcat_b": "b"})

# Creating a second data frame. Same content with dataframe 1 but no duplicates in the rxnID
data2 = data1.drop_duplicates(['rxnID', 'direction']).reset_index()

# Helpful functions
def find_average(subdf):
    mean = subdf.mean()

    return mean

def create_enz_complex_id(subdf):
    enz_complex_list = subdf.fillna('').to_list()

    enz_complex_id = "_".join(enz_complex_list)

    return enz_complex_id

rxn_names = data1['rxnID'].drop_duplicates().to_list() # List of reaction names from data1

for i in range(0,len(rxn_names)):
    mean_kcat_f = find_average(data1[(data1['rxnID'] == rxn_names[i]) & (data1['direction'] == 'f')]['kcat'])
    mean_kcat_b = find_average(data1[(data1['rxnID'] == rxn_names[i]) & (data1['direction'] == 'b')]['kcat'])
    mean_molmass = find_average(data1[data1['rxnID'] == rxn_names[i]]['molMass'])
    enz_complex_id = create_enz_complex_id(data1[data1['rxnID'] == rxn_names[i]]['uniprotID'].drop_duplicates())

    if data2['direction'].iloc[i] == 'f':
        data2.loc[i, 'kcat'] = mean_kcat_f

    else:
        data2.loc[i, 'kcat'] = mean_kcat_b

    data2.loc[i, 'molMass'] = mean_molmass
    data2.loc[i, 'uniprotID'] = enz_complex_id


data2 = data2.dropna(subset=['kcat'])
mcpam_data = data2

with pd.ExcelWriter(core_data_path, engine='openpyxl', mode='a') as writer:
    # Write the new DataFrame to a new sheet
    mcpam_data.to_excel(writer, sheet_name='mcPAM_data_core_new_structure', index=True)




