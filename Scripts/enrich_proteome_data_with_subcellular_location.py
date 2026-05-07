import pandas as pd

proteome_df = pd.read_excel('Data/E_coli_proteomics_data.xlsx', sheet_name='merged_data') # Data obtained from Shannara Kayleigh Taylor Parkins
protein_location_df = pd.read_excel('Data/uniprotkb_cell_inner_membrane_proteins_ecoli_2026_05_07.xlsx', sheet_name="raw")

# Extracting the true subcellular location from column 'Subcellular location [CC]'
protein_location_df["Cellular protein location"] = (
    protein_location_df["Subcellular location [CC]"].str.extract(r"(Cell inner membrane)")
)

# Extracting gene number with pattern b#### from 'Gene Names'
protein_location_df['Bnumber'] = (
    protein_location_df['Gene Names'].str.extract(r"(b\d{4})")
)
protein_location_df = protein_location_df[['Entry', 'Cellular protein location', 'Bnumber']]

proteome_df = pd.merge(left=proteome_df,
                       right=protein_location_df,
                       how='left',
                       left_on='uniprotID',
                       right_on='Entry')
proteome_df['Cellular protein location'] = proteome_df['Cellular protein location'].fillna(proteome_df['Type']) 
proteome_df = proteome_df.replace('cytosolic protein', 'Cytoplasm')

with pd.ExcelWriter('Data/E_coli_proteomics_data.xlsx', engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    proteome_df.to_excel(writer, sheet_name='merged_data_with_location', index=False)