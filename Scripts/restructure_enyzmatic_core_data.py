from src.PAModelpy.utils.pam_generation import _get_genes_for_proteins, merge_enzyme_complexes
import pandas as pd
import cobra
import os

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
protein2gene, gene2protein = _get_genes_for_proteins(enzyme_db, model)

new_enzyme_db = merge_enzyme_complexes(enzyme_db, gene2protein)
print(new_enzyme_db[new_enzyme_db['rxn_id'] == '14GLUCANabcpp'])