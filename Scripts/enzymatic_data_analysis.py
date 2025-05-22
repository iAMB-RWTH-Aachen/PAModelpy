import pandas as pd
from src.PAModelpy.utils.pam_generation import set_up_pam
from src.Protein import Protein

### Helper functions ###
def find_var_for_complex(enz_complex, df):
    enzymes = enz_complex.split('_')
    p_conc_in_complex = []
    alpha_in_complex = []

    for enz in enzymes:
        if df['uniprotID'].isin([enz]).any():
            p_conc_in_complex.append(df[df['uniprotID'] == enz]['SP4+TX100_avg'].values)
            alpha_in_complex.append(df[df['uniprotID'] == enz]['alpha_helix_units'].values)

    if len(p_conc_in_complex) == 0 or len(alpha_in_complex) == 0:
        p_conc_for_complex = 0
        alpha_for_complex = 0
    else:
        p_conc_for_complex = min(p_conc_in_complex)
        alpha_for_complex = max((alpha_in_complex))

    return p_conc_for_complex, alpha_for_complex

## MAIN ##

if __name__ == "__main__":
    # Generate the model
    pam_info_path = 'Results/PAM_parametrizer/Enzymatic_files/2025_05_14/proteinAllocationModel_EnzymaticData_iML1515_2.xlsx'
    pam = set_up_pam(pam_info_file=pam_info_path, sensitivity=False, membrane_sector=True)

    # Extracting the membrane enzyme complex ids from the dictionary of membrane proteins from the model
    memprot_dict = pam.sectors.get_by_id('MembraneSector').membrane_proteins
    enz_complex_list = list(memprot_dict.keys())

    # Get the enzymes from the experimental data
    data_path = "Data/E_coli_proteomics_data.xlsx"
    df = pd.read_excel(data_path, sheet_name="merged_data")

    # Initial total occupied area is 0 because no protein is allocated to the membrane yet
    total_occupied_area = 0

    for enz_complex in enz_complex_list:
        p_conc_for_complex, alpha_for_complex = find_var_for_complex(enz_complex, df)
        protein_object = Protein()
        occupied_area_per_protein = protein_object.calculate_protein_area(p_conc_for_complex, alpha_for_complex)
        total_occupied_area += occupied_area_per_protein

    print(f'Membrane occupancy [%]: {total_occupied_area / 12.10 * 100} %')