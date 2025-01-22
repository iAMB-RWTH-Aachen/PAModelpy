import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data1_path = os.path.join('Data/mcPAM_iML1515_EnzymaticData.xlsx')
data2_path = os.path.join('Data/proteinAllocationModel_mc-core_EnzymaticData_241209_multi.xlsx')
data1 = pd.read_excel(data1_path, sheet_name='mcPAM_data_core')
data2 = pd.read_excel(data2_path, sheet_name="ActiveEnzymes")

data1 = data1.drop_duplicates(subset=["rxnID"])
data1 = data1[['rxnID', 'kcat_f', 'kcat_b']]

# Unpivot the DataFrame
data1 = data1.melt(
    id_vars=['rxnID'],                # Columns to keep as is
    value_vars=['kcat_f', 'kcat_b'],  # Columns to unpivot
    var_name='direction',             # Name of the new column for direction
    value_name='kcat'                 # Name of the new column for kcat values
)

data1 = data1.dropna(subset=['kcat']).reset_index(drop=True)

# Map directions to 'f' and 'b'
data1['direction'] = data1['direction'].map({'kcat_f': 'f', 'kcat_b': 'b'})

# Drop rows where kcat is NaN
data1 = data1.dropna(subset=['kcat']).reset_index(drop=True)

merged_data = pd.merge(data2, data1, left_on=['rxn_id', 'direction'], right_on=['rxnID', 'direction'], how='right')
merged_data = merged_data.drop(columns=['rxnID', 'kcat_values']).rename(columns={'kcat': 'kcat_values'})

with pd.ExcelWriter(data2_path, engine='openpyxl', mode='a') as writer:
    # Write the new DataFrame to a new sheet
    merged_data.to_excel(writer, sheet_name='mcPAM_data_core', index=True)


