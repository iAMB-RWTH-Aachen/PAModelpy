import pandas as pd
import re

# Load the enzyme diagnostics file, adjusting column names
diagnostics_data_path = 'Results/PAM_parametrizer/Files/2025_02_20/pam_parametrizer_diagnostics_mciML1515_1.xlsx'
enzyme_db = pd.read_excel(diagnostics_data_path, sheet_name='Best_Individuals') # this sheet contains the updated kcat of each enzyme for 10 runs

# Group the enzymes after its id and the reaction id and then select only the latest run from each group
grouped = enzyme_db.groupby(['enzyme_id', 'rxn_id'])
updated_kcat_db = grouped.last().reset_index()

# Function to extract only the main reaction ID
def extract_reaction_id(reaction):
    reaction = re.sub(r'^CE_', '', reaction)   # Remove 'CE_' prefix
    reaction = re.sub(r'_[A-Z0-9]+', '', reaction)  # Remove extra identifiers
    return reaction

# Apply function to each row to extract only the main reaction ID
updated_kcat_db['rxn_id'] = updated_kcat_db['rxn_id'].apply(extract_reaction_id)

# Merge updated kcats with mcpam enzymatic data
## Load the mcpam enzymatic data
mcpam_data_path = 'Results/PAM_parametrizer/Files/2025_02_20/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
mcpam_db = pd.read_excel(mcpam_data_path, sheet_name="ActiveEnzymes")
merged_db = pd.merge(mcpam_db, updated_kcat_db, on=['enzyme_id', 'rxn_id', 'direction'], how='left')

for i, row in merged_db.iterrows():
    if not pd.isna(row["kcat[s-1]"]):
        merged_db['kcat_values'].iloc[i] = row["kcat[s-1]"]

# Drop unnecessary columns from the merged db
new_db = merged_db.drop(columns=['run_id', 'kcat[s-1]', 'ga_error'])

#Write excel datasheet
new_db_data_path = 'Results/PAM_parametrizer/Files/2025_02_20/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'
with pd.ExcelWriter(new_db_data_path, engine='openpyxl', mode='a') as writer:
    # Write the new DataFrame to a new sheet
    new_db.to_excel(writer, sheet_name='diagnostics_1', index=True)