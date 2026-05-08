import pandas as pd

proteome_df = pd.read_excel('Data/E_coli_proteomics_data.xlsx', sheet_name='merged_data_with_location') # Data obtained from Shannara Kayleigh Taylor Parkins
uniprot_protein_df = pd.read_excel('Data/uniprotkb_ecoli_whole_proteome_2024_05_02.xlsx', sheet_name="Sheet0")

# Extracting the true subcellular location from column 'Subcellular location [CC]'
uniprot_protein_df["Cellular protein location"] = (
    uniprot_protein_df["Subcellular location [CC]"].str.extract(r"(Cell inner membrane)")
)

# Extracting gene number with pattern b#### from 'Gene Names'
uniprot_protein_df['Bnumber'] = (
    uniprot_protein_df['Gene Names'].str.extract(r"(b\d{4})")
)
uniprot_protein_df = uniprot_protein_df[['Entry', 'Cellular protein location', 'Bnumber']]

proteome_df = pd.merge(left=proteome_df,
                       right=uniprot_protein_df,
                       how='left',
                       left_on='uniprotID',
                       right_on='Entry')
proteome_df['Cellular protein location'] = proteome_df['Cellular protein location'].fillna(proteome_df['Type']) 
proteome_df = proteome_df.replace('cytosolic protein', 'Cytoplasm')

with pd.ExcelWriter('Data/E_coli_proteomics_data.xlsx', engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    proteome_df.to_excel(writer, sheet_name='merged_data_with_location', index=False)