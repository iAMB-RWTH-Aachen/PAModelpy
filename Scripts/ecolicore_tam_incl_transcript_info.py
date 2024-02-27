import pandas as pd
import os

from Scripts.tam_generation import set_up_toy_tam, set_up_ecolicore_tam
sinha_ref_conditions = {
    'Holm et al': ['REF', 'NOX', 'ATP'], #WT and NADH  overexpression conditions, mu 0.72, 0.65,0.58 h-1 respectively
    'Ishii et al': ['WT_0.7h-1'],
    'Gerosa et al.': ['Glycerol','Glucose','Acetate', 'Pyruvate','Gluconate','Succinate','Galactose','Fructose'],
}
TRANSCRIPT_FILE_PATH = os.path.join('Data', 'TAModel', 'Sinha-etal_2021_transcript-data.xlsx')


mrna_vs_mu_slope = 2.64E-10
mrna_vs_mu_intercept = 3.59E-11


def get_transcript_data(transcript_file_path:str = TRANSCRIPT_FILE_PATH, mmol = True,
                        reference: str = 'Holm et al', growth_rates = [0.72, 0.65,0.58]):
    expression_data = pd.read_excel(transcript_file_path, sheet_name=reference, index_col=0)
    # Normalize columns by dividing by the sum of each column
    expression_data_normalized = expression_data.div(expression_data.sum(axis=0), axis=1)
    if mmol:
        mrna_conc = [mrna_vs_mu_intercept + mrna_vs_mu_slope*mu for mu in growth_rates]
        expression_data_mmol = expression_data_normalized.apply(lambda row: row * mrna_conc, axis=1)
        return expression_data_mmol
    else:
        return expression_data_normalized

if __name__ == '__main__':
    tam = set_up_ecolicore_tam()
    tam.optimize()
    transcript_data_mmol = get_transcript_data()

    for gene, expression_data in transcript_data_mmol.iterrows():
        transcript_id = 'mRNA_'+gene
        if not transcript_id in tam.transcripts: continue
        transcript = tam.transcripts.get_by_id('mRNA_'+gene)
        #testing wildtype condition
        transcript.change_concentration(concentration=expression_data[0],
                                        error= expression_data[0]*0.01)
    print(tam.reactions.get_by_id('EX_glc__D_e'))

    tam.optimize()
    print(tam.summary())
