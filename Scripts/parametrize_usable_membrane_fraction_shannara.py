import pandas as pd
from src.PAModelpy.utils.pam_generation import set_up_pam, set_up_core_pam
from src.Protein import Protein

### Helper functions ###
def find_var_for_complex(enz_complex, df):
    enzymes = enz_complex.id.split('_')
    p_conc_in_complex = []
    alpha_in_complex = []

    for enz in enzymes:
        row = df[df['uniprotID'] == enz]
        if not row.empty and row["Cellular protein location"].iloc[0] == "Cell inner membrane":
            p_conc_in_complex.append(df[df['uniprotID'] == enz]['SP4+TX100_avg'].values)
            alpha_in_complex.append(df[df['uniprotID'] == enz]['alpha_helix_units'].values)

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
    data_path = "Data/proteome_data_full_shannara2024.xlsx"
    df = pd.read_excel(data_path, sheet_name="merged_data_with_location")
    count_memprot_not_inner = df[(df["Cellular protein location"] == "membrane protein")].shape[0]
    count_inner_memprot = df[(df["Cellular protein location"] == "Cell inner membrane")].shape[0]
    count_cyto_prot = df[(df["Cellular protein location"] == "Cytoplasm")].shape[0]
    print(f'Count of membrane proteins that are not inner membrane protein: {count_memprot_not_inner}')
    print(f'Count of membrane proteins that are inner membrane protein: {count_inner_memprot}')
    print(f'Count of cytosolic proteins: {count_cyto_prot}')

    # Initial total occupied area is 0 because no protein is allocated to the membrane yet
    total_occupied_area = 0

    for enz_complex in enzyme_complex_list:
        p_conc_for_complex, alpha_for_complex = find_var_for_complex(enz_complex, df)
        protein_object = Protein()
        protein_object.alpha_area = 1.4 * 1e-18 # m2 
        occupied_area_per_protein = protein_object.calculate_protein_area(p_conc_for_complex, alpha_for_complex)
        total_occupied_area += occupied_area_per_protein

    print(f'Membrane occupancy [%]: {total_occupied_area / 10.68 * 100} %') # 10.68 um2 is the inner membrane area for ecoli k-12 mg1655 at growth rate 0.67, grown in minimal glucose media