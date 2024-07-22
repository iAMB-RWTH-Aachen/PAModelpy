import os
import pandas as pd

PAM_DATA_FILE_PATH = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls')
RESULT_DATA_FILE = os.path.join('Results', 'protein_costs_glycolysis_tca.xlsx')

active_enzyme_info = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes')

tca_glycolysis_reactions = ['GLCptspp','PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PDH','ATPS4rpp',
                            'CYTBO34pp','NADH16pp','CS','ACONTa','ACONTb','ICDHyr','AKGDH','SUCOAS','FUM','MDH',
                            'SUCDi','NADTRHD','ACKr', 'ACt2rpp','PTAr', 'ACACT1r','ACACT2r','ACACT3r','ACACT4r',
                            'ACACT7r','ACACT8r']

protein_efficiency_df = active_enzyme_info.filter(['rxnID','molMass','kcat'])
protein_efficiency_df = protein_efficiency_df.assign(rxnName = lambda x: x.rxnID.str.rstrip('_fb'))
protein_efficiency_glyc_tca = protein_efficiency_df[protein_efficiency_df.rxnName.isin(tca_glycolysis_reactions)]
protein_efficiency_glyc_tca  = protein_efficiency_glyc_tca.assign(kcat = lambda x:x['kcat']*3600,
                                                                  relEfficiency = lambda x: x['molMass']/x['kcat']
                                                                  ).sort_values(by = 'relEfficiency', ascending=False)
print(protein_efficiency_glyc_tca)

protein_efficiency_glyc_tca.to_excel(RESULT_DATA_FILE, sheet_name='protein_costs')