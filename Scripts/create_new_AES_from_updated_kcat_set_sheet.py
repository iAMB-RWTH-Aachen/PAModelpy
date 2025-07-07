import pandas as pd
import os 
import numpy as np

from src.PAModelpy import CatalyticEvent

UPDATED_KCAT_FILE = os.path.join('Results', 'From_kcat_dataset_20250627', 'protein_occupancy_data_updated_kcats.xlsx')
NEW_AES_SUFFIX = 'iML1515_20250702'
ENZYMATIC_FILE_OLD = os.path.join('Data', 'mcPAM_iML1515_EnzymaticData_250627_TE_UE_modified.xlsx')

def reshape_flux_directions(df):

    # Melt the dataframe
    df_melted = pd.melt(df,
                        id_vars=[col for col in df.columns if col not in ['Forward Flux', 'Backward Flux']],
                        value_vars=['Forward Flux', 'Backward Flux'],
                        var_name='direction',
                        value_name='kcat_values')

    # Remove ' Flux' and lowercase
    df_melted['direction'] = df_melted['direction'].str.replace(' Flux', '').str.lower()

    # Replace full words with abbreviations
    df_melted['direction'] = df_melted['direction'].replace({
        'forward': 'f',
        'backward': 'b'
    })

    return df_melted

def search_index_in_parameter_file(df:pd.DataFrame, protein:str, reaction:str, direction:str):
    all_protein_rows = df[df.enzyme_id == protein]
    all_rxn_and_protein_rows = all_protein_rows[all_protein_rows.rxn_id == reaction]
    the_row = all_rxn_and_protein_rows[all_rxn_and_protein_rows.direction == direction]
    return the_row.index

def create_new_aes_parameter_file(old_enzymatic_file,
                                  updated_kcat_df:pd.DataFrame,
                                  new_aes_suffix: str = NEW_AES_SUFFIX):

    # extract reaction id from catalytic event id
    updated_kcat_df['Reaction'] = [CatalyticEvent._extract_reaction_id_from_catalytic_reaction_id(id) for id in
                                         updated_kcat_df['Reaction']]
    old_aes_df = old_enzymatic_file['ActiveEnzymes']

    for index, row in updated_kcat_df.iterrows():
        aes_index = search_index_in_parameter_file(old_aes_df, row['enzyme_id'], row['Reaction'], row['direction'])
        old_aes_df.loc[aes_index, 'kcat_values'] = row['kcat_values']

    write_mode = 'w'
    kwargs = {}
    result_enzymatic_file = os.path.join('Results', 'From_kcat_dataset_20250627', 'result_enzymatic_files',
                               f'proteinAllocationModel_EnzymaticData_{new_aes_suffix}.xlsx')
    if os.path.isfile(result_enzymatic_file):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(result_enzymatic_file,
            mode=write_mode, engine='openpyxl', **kwargs) as writer:
        old_aes_df.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
        for sheet, df in old_enzymatic_file.items():
            if sheet != 'ActiveEnzymes':
                df.to_excel(writer, sheet_name = sheet, index=False)

if __name__ == '__main__':
    # Load the dfs
    updated_kcat_df = pd.read_excel(UPDATED_KCAT_FILE, sheet_name='enzymatic_file_250627')
    old_enzymatic_file = pd.read_excel(ENZYMATIC_FILE_OLD, sheet_name=None)
    
    # Reshape the columns "Foward Flux" and "Backward Flux"
    reshaped_updated_kcat_df = reshape_flux_directions(updated_kcat_df)

    create_new_aes_parameter_file(old_enzymatic_file=old_enzymatic_file,
                                  updated_kcat_df=reshaped_updated_kcat_df,
                                  new_aes_suffix=NEW_AES_SUFFIX)
    

    
