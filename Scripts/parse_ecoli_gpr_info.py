from cobra.io import read_sbml_model
from cobra import Reaction
import os
import pandas as pd
import numpy as np
from datetime import date

TAM_DATA_PATH = os.path.join('Data', 'TAModel')
TAM_DATA_FILE = os.path.join(TAM_DATA_PATH, str(date.today()) +'_gene_enzyme_reaction_relation_Ecoli.xlsx')

def parse_gpr_relationships_from_Ecocyc():
    # files last downloaded from Ecocyc: 2024 - 01 - 25
    # manual curation of the tabs in the gene file is required!
    enzymes_to_genes = pd.read_csv(os.path.join(TAM_DATA_PATH , 'All-enzymes-of-E.-coli-K-12-substr.-MG1655.tsv'), sep = '\t')
    genes_info = pd.read_csv(os.path.join(TAM_DATA_PATH, 'All-genes-of-E.-coli-K-12-substr.-MG1655.tsv'), sep = '\t').drop('Object ID', axis = 1)
    genes_info['Accession-1'] = genes_info['Accession-1'].str.strip()

    enzyme_gene_info = pd.merge(enzymes_to_genes, genes_info, how='outer', left_on = 'Gene', right_on= 'Gene Name')
    enzyme_gene_info['mrna_length'] = enzyme_gene_info['Right-End-Position']-enzyme_gene_info['Left-End-Position']

    enzyme_gene_info = enzyme_gene_info.drop(['Catalyzes', 'Sequence - polypeptide sequence', 'Gene Name',
                                              'Left-End-Position', 'Right-End-Position', 'Product'], axis=1)
    model_path = os.path.join('Models', 'iML1515.xml')
    model = read_sbml_model(model_path)


    enzyme_gene_reaction_relation = pd.DataFrame(columns = list(enzyme_gene_info.columns)+ ['Reaction', 'gpr'])

    #match gpr relationships and average the molecular mass
    for index, row in enzyme_gene_info.iterrows():
        molmass = row['Molecular-Weight-KiloDaltons']
        if not isinstance(molmass, float):
            row['Molecular-Weight-KiloDaltons'] = np.mean([float(mass) for mass in molmass.split('//')])

        if row['Accession-1'] in model.genes:
            gene = model.genes.get_by_id(row['Accession-1'])
            for rxn in gene.reactions:
                gpr_info = parse_and_or_gpr_relations_from_string(rxn.gpr.to_string())

                enzyme_gene_reaction_relation.loc[len(enzyme_gene_reaction_relation)] = row.to_list() + [rxn.id] + [gpr_info]

    #find model reactions which are not matched and not exchange reactions
    not_matched_reactions = []
    for rxn in model.reactions:
        if rxn.id not in enzyme_gene_reaction_relation['Reaction'].to_list() and 'EX' not in rxn.id:
            not_matched_reactions.append(rxn.id)

    #make df look pretty
    enzyme_gene_reaction_relation = enzyme_gene_reaction_relation[
        ['Gene', 'Enzyme', 'Reaction', 'Accession-1', 'Object ID',
         'Molecular-Weight-KiloDaltons', 'mrna_length', 'gpr']]
    enzyme_gene_reaction_relation.columns = ['Gene', 'Enzyme', 'Reaction', 'gene_id', 'enzyme_id',
                                             'molmass_kDa', 'mrna_length', 'gpr']

    #write to excel
    with pd.ExcelWriter(TAM_DATA_FILE) as writer:
        enzyme_gene_reaction_relation.to_excel(writer, sheet_name='enzyme-gene-reaction')
        pd.DataFrame({'not_matched': not_matched_reactions}).to_excel(writer, sheet_name='non-matched_reactions')


def parse_and_or_gpr_relations_from_string(gpr_info:str):
    gpr_list = []
    for gene_relation in gpr_info.split(') or ('):
        gpr_list.append(gene_relation.split(' and '))
    return gpr_list

def parse_enzymatic_data_information():
    enzyme_gene_reaction_relation = pd.read_excel(TAM_DATA_FILE, sheet_name='enzyme-gene-reaction')
    enzyme_info = pd.read_excel(os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls'),
                                sheet_name='ActiveEnzymes')

    tam_info = pd.merge(enzyme_info, enzyme_gene_reaction_relation, how = 'left', left_on='rxnID', right_on='Reaction')




if __name__ == '__main__':
    parse_gpr_relationships_from_Ecocyc()