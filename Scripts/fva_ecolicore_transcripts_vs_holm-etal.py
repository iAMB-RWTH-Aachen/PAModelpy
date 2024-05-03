import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import cobra

from src.flux_analysis import flux_variability_analysis
from src.TAModelpy import Transcript

from Scripts.tam_generation import set_up_toy_tam, set_up_ecolicore_tam
from Scripts.ecolicore_tam_incl_transcript_info import get_transcript_data

TRANSCRIPT_FILE_PATH = os.path.join('Data', 'TAModel', 'Sinha-etal_2021_transcript-data.xlsx')


ecolicore_tam = set_up_ecolicore_tam()
# result = flux_variability_analysis(ecolicore_tam, variable_type = Transcript, variable_list = [trans for trans in ecolicore_tam.transcripts])
# result.to_excel('Results/fva_transcripts_higher_ub.xlsx')
result = pd.read_excel('Results/fva_transcripts_higher_ub.xlsx',index_col=0)#.reset_index(inplace=True)
result = result.reset_index()
result['index'] = result['index'].apply(lambda x: x.split('_')[1:])

transcript_data = get_transcript_data()

result_dict = {}
for index, row in result.iterrows():
    genes = row['index']
    min_val = row['minimum']
    max_val = row['maximum']
    for gene in genes:
        transcript_values = []
        print(gene)
        if gene in transcript_data.index:
            transcript_values.append(transcript_data.loc[gene].REF)
    if len(transcript_values)>0:
        ref_val = np.min(transcript_values)
        result_dict['_'.join(genes)] = [min_val <= ref_val*1e-6 <= max_val, min_val, ref_val, max_val]

result_df = pd.DataFrame.from_dict(result_dict, orient='index', columns=['isinrange', 'min', 'reference', 'max'])


print(result_df.to_markdown())
#TODO have to further validate this and visualize the results

