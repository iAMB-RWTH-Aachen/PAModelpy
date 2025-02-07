from src.PAModelpy.utils.pam_generation import get_protein_gene_mapping, merge_enzyme_complexes
import pandas as pd
import cobra
import os

# Helpful functions
def find_average(subdf):
    mean = subdf.mean()

    return mean

def sum_value(subdf):
    total = subdf.sum()

    return total

# Load the enzyme database, adjusting column names
data_path = 'Data/mcPAM_iML1515_EnzymaticData.xlsx'
enzyme_db = pd.read_excel(data_path, sheet_name='mcPAM_data_core')
enzyme_db = enzyme_db.rename(columns={'rxnID': 'rxn_id',
                                      'uniprotID': 'enzyme_id',
                                      'm_gene': 'gene',
                                      'm_gene_reaction_rule': 'GPR'})

# Load the model
model = cobra.io.read_sbml_model(os.path.join('Models', 'iML1515.xml'))

# Create rxn2protein and protein2gpr dictionary
protein2gene, gene2protein = get_protein_gene_mapping(enzyme_db, model)

new_enzyme_db = merge_enzyme_complexes(enzyme_db, gene2protein)

# Unpivot the DataFrame 'kcat_f/b' --> 'kcat' and 'direction (f/b)'
columns = ['enzyme_id', 'rxn_id', 'molMass', 'GPR', 'gene', 'direction', 'kcat']
new_enzyme_db = new_enzyme_db.melt(
    id_vars=columns[0:5],                # Columns to keep as is
    value_vars=['kcat_f', 'kcat_b'],  # Columns to unpivot
    var_name='direction',             # Name of the new column for direction
    value_name='kcat'                 # Name of the new column for kcat values
)
new_enzyme_db["direction"] = new_enzyme_db["direction"].replace({"kcat_f": "f", "kcat_b": "b"})

# Remove duplicates in rxn_ID and average out the alpha number
new_enzyme_db_copy = new_enzyme_db.drop_duplicates(['rxn_id', 'direction']).reset_index()

rxn_names = new_enzyme_db['rxn_id'].drop_duplicates().to_list() # List of reaction names from data1

for i in range(0,len(rxn_names)):
    mean_kcat_f = find_average(new_enzyme_db[(new_enzyme_db['rxn_id'] == rxn_names[i]) & (new_enzyme_db['direction'] == 'f')]['kcat'])
    mean_kcat_b = find_average(new_enzyme_db[(new_enzyme_db['rxn_id'] == rxn_names[i]) & (new_enzyme_db['direction'] == 'b')]['kcat'])
    mean_molmass = find_average(new_enzyme_db[new_enzyme_db['rxn_id'] == rxn_names[i]]['molMass'])

    if new_enzyme_db_copy['direction'].iloc[i] == 'f':
        new_enzyme_db_copy.loc[i, 'kcat'] = mean_kcat_f

    else:
        new_enzyme_db_copy.loc[i, 'kcat'] = mean_kcat_b

    new_enzyme_db_copy.loc[i, 'molMass'] = mean_molmass

new_enzyme_db_copy = new_enzyme_db_copy.dropna(subset=['kcat'])
mcpam_data = new_enzyme_db_copy

# Write excel datasheet
core_data_path = os.path.join('Data/proteinAllocationModel_mc-core_EnzymaticData_241209_multi.xlsx')
with pd.ExcelWriter(core_data_path, engine='openpyxl', mode='a') as writer:
    # Write the new DataFrame to a new sheet
    mcpam_data.to_excel(writer, sheet_name='mcPAM_data_core_new_structure_w', index=True)