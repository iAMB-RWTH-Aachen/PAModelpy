import pandas as pd
from src.PAModelpy.utils.pam_generation import set_up_pam
from src.Protein import Protein

### Helper functions ###
def find_var_for_complex(enz_complex, df):
    enzymes = enz_complex.id.split('_')
    p_conc_in_complex = []
    alpha_in_complex = []

    for enz in enzymes:
        row = df[df['enzyme_id'] == enz]
        if not row.empty and row["Cellular protein location"].iloc[0] == "Cell inner membrane":
            p_conc_in_complex.append(df[df['enzyme_id'] == enz].max(axis=1).mean())
            alpha_in_complex.append(df[df['enzyme_id'] == enz]['alpha_numbers'].mean().mean())

    if len(p_conc_in_complex) == 0 or len(alpha_in_complex) == 0:
        p_conc_for_complex = 0
        alpha_for_complex = 0
    else:
        p_conc_for_complex = min(p_conc_in_complex)
        alpha_for_complex = sum(alpha_in_complex)

    return p_conc_for_complex, alpha_for_complex

## MAIN ##

if __name__ == "__main__":
    # Generate the model
    pam_info_path = 'Data/proteinAllocationModel_EnzymaticData_iML1515_10.xlsx'
    pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False)
    pam.optimize()
    print(pam.objective.value)

    # Extracting enzymes from the model
    enzyme_complex_list = pam.enzyme_variables

    # Get the enzymes from the experimental data
    proteome_df_extract = pd.read_excel('Data/proteome_data_extract_schmidt2016.xlsx',
                            sheet_name='ProteinCopies',
                            engine='openpyxl')

    proteome_df_full = pd.read_excel('Data/proteome_data_full_schmidt2016.xlsx',
                                sheet_name='Table S6',
                                engine='openpyxl',
                                header=2)

    protein2location_df = pd.read_excel('Data/proteome_data_extract_schmidt2016.xlsx',
                                sheet_name='SubcellularLocation',
                                engine='openpyxl')
    
    membrane_enzymes_df = pd.read_excel('Data/proteinAllocationModel_EnzymaticData_iML1515_10.xlsx',
                                sheet_name='MembraneEnzymes', usecols=['enzyme_id', 'alpha_numbers'],
                                engine='openpyxl').drop_duplicates()
    
    proteome_df_full = proteome_df_full[['Bnumber', 'Uniprot Accession', 'Glucose']]
    # proteome_df_full = proteome_df_full.rename(columns = {'Glucose.2': 'Glucose'})

    proteome_df = pd.merge(left=proteome_df_extract,
                        right=proteome_df_full,
                        how='left',
                        on=['Bnumber', 'Glucose']) # in fg/cell
    proteome_df = pd.merge(left=proteome_df,
                        right=protein2location_df,
                        how='left',
                        on="Uniprot Accession")
    proteome_df = pd.merge(left=proteome_df,
                        right=membrane_enzymes_df,
                        how='left',
                        left_on="Uniprot Accession",
                        right_on='enzyme_id')

    # Move this column to the front
    col = 'Uniprot Accession'
    proteome_df = proteome_df[[col] + [c for c in proteome_df.columns if c != col]]
    proteome_df = proteome_df.set_index(["Bnumber", "Uniprot Accession"])

    # Initial total occupied area is 0 because no protein is allocated to the membrane yet
    total_occupied_area = 0

    for enz_complex in enzyme_complex_list:
        p_conc_for_complex, alpha_for_complex = find_var_for_complex(enz_complex, proteome_df)
        protein_object = Protein()
        protein_object.alpha_area = 1.4 * 1e-18 # m2 
        occupied_area_per_protein = protein_object.calculate_protein_area(p_conc_for_complex, alpha_for_complex, 'protein copies/cell')
        total_occupied_area += occupied_area_per_protein

    print(f'Membrane occupancy [%]: {total_occupied_area / 10.68 * 100} %') # 10.68 um2 is the inner membrane area for ecoli k-12 mg1655 at growth rate 0.67, grown in minimal glucose media