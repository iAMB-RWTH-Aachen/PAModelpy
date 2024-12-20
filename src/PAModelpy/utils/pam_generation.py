import pandas as pd
import numpy as np
import cobra
from typing import TypedDict, Literal, Union, Tuple, Iterable
import re
import os

from collections import defaultdict
from dataclasses import dataclass, field


from ..PAModel import PAModel
from ..EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from ..MembraneSector import MembraneSector
from ..configuration import Config

DEFAULT_MOLMASS = 39959.4825 #kDa
DEFAULT_KCAT = 11 #s-1

class EnzymeInformation(TypedDict):
    enzyme_id:str
    f:float  # Forward kcat
    b:float  # Backward kcat
    genes: list[str]
    protein_reaction_association: str
    molmass: float

@dataclass
class ReactionInformation:
    rxnid: str
    enzymes: defaultdict[str,EnzymeInformation]= field(default_factory=defaultdict)
    model_reactions: Iterable[cobra.Reaction] = None

    def get_reaction_from_model(self, model):
        if self.rxnid in model.reactions:
            self.model_reactions = [model.reactions.get_by_id(self.rxnid)]
        else:
            self.model_reactions = [rxn for rxn in model.reactions.query(self.rxnid) if '_copy' in rxn.id]
        return self.model_reactions

    def check_reaction_reversibility(self) -> Literal[-1,0,1]:
        if self.model_reactions is None:
            raise AttributeError('Please first set the model reaction using ReactionInformation.get_reaction_from_model(cobra.Model)')
        if all([rxn.lower_bound>=0 for rxn in self.model_reactions]):
            # irreversible reaction (forward direction)
            rev = 1
        elif all([rxn.upper_bound<=0 for rxn in self.model_reactions]):
            # irreversible reaction (reverse direction)
            rev = -1
        else:
            # reversible reaction
            rev = 0
        return rev

    def to_dict(self):
        return {self.rxnid: dict(self.enzymes)}

def enzyme_information(rxn_id: str,
                       genes: list[str],
                       protein_reaction_association: list[str] = None,
                       enzyme_id: str = None,
                       kcat_f: float= DEFAULT_KCAT,
                       kcat_b: float= DEFAULT_KCAT,
                       molmass:float = DEFAULT_MOLMASS) -> EnzymeInformation:
    if enzyme_id is None:
        enzyme_id = 'Enzyme_' + rxn_id
        protein_reaction_association = [[enzyme_id]]
    if molmass is None or np.isnan(molmass): molmass = DEFAULT_MOLMASS

    return {'enzyme_id': enzyme_id,
            'f':_check_kcat_value(kcat_f),
            'b':_check_kcat_value(kcat_b),
            'genes':genes,
            'protein_reaction_association':protein_reaction_association,
            'molmass':molmass
            }

def parse_gpr_information(gpr_info:str,
                          genes:list=None,
                          enzyme_id:str=None,
                          gene2protein: dict[str, str] = None) -> tuple[list,list]:
    #filter out nan entries
    if not isinstance(gpr_info, str):
        return None, None

    # #only get the genes associated with this enzyme
    gpr_list = _parse_gpr(gpr_info)
    gpr_list = _filter_sublists(gpr_list, genes)

    if genes is None: return gpr_list

    #convert the genes to the associated proteins
    enzyme_relations = []
    if '_'in enzyme_id:
        enzyme_relations = [enzyme_id.split('_')]
    for sublist in gpr_list:
        enz_sublist = []
        for item in sublist:
            if item in gene2protein.keys():
                if '_' not in gene2protein[item]:
                    enz_sublist.append(gene2protein[item])
                    enzyme_relations += [enz_sublist]
                elif gene2protein[item].split('_') not in enzyme_relations:
                    enzyme_relations += [gene2protein[item].split('_')]
        enzyme_relations = _filter_sublists(enzyme_relations, enzyme_id.split('_'))
    return sorted(gpr_list), sorted(enzyme_relations)

def get_protein_gene_mapping(enzyme_db: pd.DataFrame, model) -> tuple[dict, dict]:
    protein2gene = defaultdict(list)
    gene2protein = {}
    for index, row in enzyme_db.iterrows():
        # Parse data from the row
        rxn_id = row['rxn_id']
        if rxn_id not in model.reactions: continue
        rxn = model.reactions.get_by_id(rxn_id)
        # get the identifiers and replace nan values by dummy placeholders
        enzyme_id = row['enzyme_id']
        gene_id = row['gene']

        # check if there are genes associates with the reaction
        if len(rxn.genes) > 0 or isinstance(gene_id, str):
            if not isinstance(enzyme_id, str):
                enzyme_id = 'Enzyme_' + rxn_id
            if not isinstance(gene_id, str):
                gene_id = [gene.id for gene in rxn.genes][0]

            gene2protein[gene_id] = enzyme_id

            # Create gene-protein-reaction associations
            if enzyme_id not in protein2gene:
                protein2gene[enzyme_id].append(gene_id)

    return dict(protein2gene), gene2protein

def _parse_gpr(gpr_info):
    # Split the string by 'or' first
    or_groups = gpr_info.split(' or ')
    parsed_gpr = []
    for group in or_groups:
        # Within each 'or' group, split by 'and'
        and_genes = group.split(' and ')
        # Clean up each gene string
        and_genes = [gene.strip(' ()') for gene in and_genes]
        parsed_gpr.append(and_genes)
    return parsed_gpr

def _filter_sublists(nested_list:list[list[str]], target_list:list[str]) -> list[str]:
    """
    Keeps only sublists from `nested_list` that contain at least one identifier present in `target_list`.

    Args:
        nested_list (list of lists): The main list of sublists to filter.
        target_list (list of str): The flat list of identifiers to include.

    Returns:
        list of lists: The filtered nested list.
    """
    filtered_list = [sublist for sublist in nested_list if any([item in target_list for item in sublist])]

    return filtered_list



def _check_kcat_value(kcat:float):
    return [kcat if not np.isnan(kcat) else 0][0]

def _check_enzyme_id(gene_id:str,
                     enzyme_id:Union[str, float],
                     gene2protein: dict[str, str]) -> str:
    if not isinstance(enzyme_id, str):
        enzyme_id = gene2protein[gene_id]
    return enzyme_id

def _check_rxn_identifier_format(rxn_id:str) -> str:
    for substring in ['ppi', 'pp', 'ex']:
        if substring in rxn_id:
            return rxn_id.replace(substring, '')
    return rxn_id


def _extract_reaction_id(input_str):
    # Regular expression to remove 'copy' and the following digits (if present)
    match = re.match(r'([^\_]+)_copy\d*', input_str)

    if match:
        return match.group(1)
    else:
        return input_str  # If 'copy' is not found, return the original string

def _check_if_all_model_reactions_are_in_rxn_info2protein(model: cobra.Model,
                                                          rxn_info2protein:dict[str, ReactionInformation],
                                                          protein2gpr:defaultdict[str,list]
                                                          ) -> Tuple[dict[str, ReactionInformation],defaultdict[str,list]]:
    for rxn in model.reactions:
        rxn_id = _extract_reaction_id(
            rxn.id)  # some reactions ids are associated with copy numbers, only filter for the actual reaction id
        if not (
                rxn_id not in rxn_info2protein.keys()
                and 'EX'.lower() not in rxn.id.lower()  # is the reaction an exchange with the environment?
                and 'BIOMASS' not in rxn.id  # is the reaction a pseudoreaction?
                and len(rxn._genes) > 0  # is the reaction associated with enzymes?
                and list(rxn._genes)[0].id != 's0001'
        ): continue

        print('No enzyme information found for reaction: ' + rxn.id)

        rxn_info = ReactionInformation(rxn.id)
        rxn_info.model_reactions = rxn
        kwargs = {}
        if rxn_info.check_reaction_reversibility() > 0:
            kwargs = {'kcat_b':0}
        elif rxn_info.check_reaction_reversibility() < 0:
            kwargs = {'kcat_f': 0}


        gpr_info = parse_gpr_information(rxn.gpr)
        if gpr_info is None: gpr_info = [[f'gene_{rxn.id}']]
        enzyme_info = enzyme_information(rxn.id, rxn._genes, **kwargs)
        rxn_info.enzymes[enzyme_info['enzyme_id']] = enzyme_info

        rxn_info2protein[rxn.id] = rxn_info
        protein2gpr[enzyme_info['enzyme_id']] = gpr_info

    return rxn_info2protein, protein2gpr

def _get_genes_for_proteins(enzyme_db: pd.DataFrame, model) -> dict:
    protein2gene = defaultdict(list)
    gene2protein = {}
    for index, row in enzyme_db.iterrows():
        # Parse data from the row
        rxn_id = row['rxn_id']
        if rxn_id not in model.reactions: continue
        rxn = model.reactions.get_by_id(rxn_id)
        # get the identifiers and replace nan values by dummy placeholders
        enzyme_id = row['enzyme_id']
        genes = row['gene']

        # check if there are genes associates with the reaction
        if (not len(rxn.genes) > 0) and (not isinstance(genes, list)): continue
        if not isinstance(enzyme_id, str):
            enzyme_id = 'Enzyme_' + rxn_id
        if not isinstance(genes, list):
            genes = [gene.id for gene in rxn.genes]

        for gene in genes:
            gene2protein[gene] = enzyme_id
        # Create gene-protein-reaction associations
        protein2gene[enzyme_id].append(genes)
        # assume that there is a single protein if the previous relation was not assigned a gene
        if 'gene_' in protein2gene[enzyme_id][0]:
             protein2gene[enzyme_id] = [genes]

    return dict(protein2gene), gene2protein

def _get_kcat_value_from_catalytic_reaction_df(catalytic_reaction_df: pd.DataFrame,
                                               direction:str) -> float:
    if direction not in catalytic_reaction_df.direction.values:
        return 0
    else:
        return catalytic_reaction_df.kcat_values[
            catalytic_reaction_df.direction == direction].iloc[0]


def parse_reaction2protein(enzyme_db: pd.DataFrame, model: cobra.Model) -> dict:
    rxn_info2protein = {}
    protein2gpr = defaultdict(list)
    #remove copy number substrings from the reaction to make it matchable to enzyme information
    filtered_model_reactions = [_extract_reaction_id(r.id) for r in model.reactions]

    # replace NaN values with unique identifiers
    enzyme_db.loc[enzyme_db['enzyme_id'].isnull(), 'enzyme_id'] = [f'E{i}' for i in
                                                                   range(enzyme_db['enzyme_id'].isnull().sum())]

    enzyme_db.loc[enzyme_db['gene'].isnull(), 'gene'] = [[f'gene_{i}'] for i in
                                                                   range(enzyme_db['gene'].isnull().sum())]

    protein2gene, gene2protein = _get_genes_for_proteins(enzyme_db, model)

    # parse the information for all gene-protein-reaction relations in the dataframe
    for (rxn_id, enzyme_id, gpr), catalytic_reaction_info in enzyme_db.groupby(
            ['rxn_id', 'enzyme_id', 'GPR']):
        # only parse those reactions which are in the model
        if rxn_id not in filtered_model_reactions: continue

        rxn_info = rxn_info2protein.setdefault(rxn_id, ReactionInformation(rxn_id))
        #sometimes, multiple copies are associated with a single reaction
        rxns = rxn_info.get_reaction_from_model(model)
        genes = catalytic_reaction_info.gene.iloc[0]
        for rxn in rxns:
            # If no genes are associated with the reaction, this reaction is not catalyzed by an enzyme
            if (not len(rxn.genes) > 0) and (not isinstance(genes, list)):  continue
            # enzyme_id = check_enzyme_id(gene_id, enzyme_id, gene2protein)

            # get the gene-protein-reaction-associations for this specific enzyme
            gene_reaction_relation, protein_reaction_relation = parse_gpr_information(gpr,
                                                                                      genes,
                                                                                      enzyme_id,
                                                                                      gene2protein)

            protein2gpr[enzyme_id]+= gene_reaction_relation

            enzyme_info = enzyme_information(rxn_id=rxn.id,
                                             genes=genes,
                                             protein_reaction_association=protein_reaction_relation,
                                             enzyme_id=enzyme_id,
                                             kcat_f=_get_kcat_value_from_catalytic_reaction_df(
                                                        catalytic_reaction_info, 'f'),
                                             kcat_b=_get_kcat_value_from_catalytic_reaction_df(
                                                        catalytic_reaction_info, 'b'),
                                             molmass=catalytic_reaction_info.molMass.iloc[0])
            rxn_info.enzymes[enzyme_id] = enzyme_info
            rxn_info2protein[rxn.id] = rxn_info

    # if no enzyme info is found, add dummy enzyme with median kcat and molmass
    rxn_info2protein, protein2gpr = _check_if_all_model_reactions_are_in_rxn_info2protein(model, rxn_info2protein, protein2gpr)

    #convert the dataobject dict to a normal dict for correct parsing by PAModelpy
    rxn2protein = {rxn_id: dict(rxn_info.enzymes) for rxn_id, rxn_info in rxn_info2protein.items()}

    return rxn2protein, dict(protein2gpr)


# Function to parse GPR and determine multimer
def merge_enzyme_complexes(df, gene2protein):
    collapsed_rows = []

    for rxn_id, group in df.groupby('rxn_id'):
        for _, row in group.iterrows():
            # Parse GPR
            gpr_list, enzyme_relations = parse_gpr_information(
                row['GPR'], row['gene'], row['enzyme_id'], gene2protein
            )
            # Collapse "and" relationships into multimer ID
            if enzyme_relations:
                for gene_list, enzyme_list in zip(gpr_list, enzyme_relations):
                    row_copy = row.copy()
                    row_copy['enzyme_id'] = "_".join(enzyme_list)  # Replace gene with multimer
                    row_copy['gene'] = gene_list  # add all the annotations to the corresponding gene

                    collapsed_rows.append(row_copy)
            else:
                collapsed_rows.append(row)

    # Create a new dataframe with collapsed rows
    collapsed_df = pd.DataFrame(collapsed_rows)
    return collapsed_df

def set_up_pam(pam_info_file:str = '',
               model:Union[str, cobra.Model] = 'Models/iML1515.xml',
               config:Config = None,
               total_protein: Union[bool, float] = True,
               active_enzymes: bool = True,
               translational_enzymes: bool = True,
               unused_enzymes: bool = True,
               membrane_sector: bool = False,
               max_membrane_area:float = 0.03,
               sensitivity:bool = True,
               enzyme_db:pd.DataFrame = None,
               adjust_reaction_ids:bool = False) -> PAModel:


    if config is None:
        config = Config()
        config.reset()

    # some other constants
    TOTAL_PROTEIN_CONCENTRATION = 0.258  # [g_prot/g_cdw]

    #setup model if a path is provided
    if isinstance(model, str):
        model = cobra.io.read_sbml_model(model)

    #check if a different total protein concentration is given
    if isinstance(total_protein, float):
        TOTAL_PROTEIN_CONCENTRATION = total_protein

    # load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        if enzyme_db is None:
            enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')
            #for some models, the reaction ids should not include 'pp' or 'ex'
            if adjust_reaction_ids:
                enzyme_db['rxn_id'] = enzyme_db['rxn_id'].apply(_check_rxn_identifier_format)
        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, model)
        active_enzyme_info = ActiveEnzymeSector(rxn2protein=rxn2protein, protein2gene=protein2gene,
                                                  configuration=config)
    else:
        active_enzyme_info = None

    if translational_enzymes:
        translational_info = pd.read_excel(pam_info_file, sheet_name='Translational')
        translation_enzyme_info = TransEnzymeSector(
            id_list=[translational_info[translational_info.Parameter == 'id_list'].loc[0, 'Value']],
            tps_0=[translational_info[translational_info.Parameter == 'tps_0'].loc[1, 'Value']],
            tps_mu=[translational_info[translational_info.Parameter == 'tps_mu'].loc[2, 'Value']],
            mol_mass=[translational_info[translational_info.Parameter == 'mol_mass'].loc[3, 'Value']],
            configuration = config)
    else:
        translation_enzyme_info = None

    if unused_enzymes:
        unused_protein_info = pd.read_excel(pam_info_file, sheet_name='UnusedEnzyme').set_index('Parameter')

        ups_0 = unused_protein_info.at['ups_0', 'Value']
        ups_mu = unused_protein_info.at['ups_mu', 'Value']

        unused_enzyme_info = UnusedEnzymeSector(
            id_list=[unused_protein_info.at['id_list', 'Value']],
            ups_mu=[ups_mu],
            ups_0=[ups_0],
            mol_mass=[unused_protein_info.at['mol_mass', 'Value']],
            configuration = config)
    else:
        unused_enzyme_info = None

    if membrane_sector:
        membrane_info = pd.read_excel(pam_info_file, sheet_name='Membrane').set_index('Parameter')
        active_membrane_info = pd.read_excel(pam_info_file, sheet_name='MembraneEnzymes').set_index('enzyme_id')

        area_avail_0 = membrane_info.at['area_avail_0','Value']
        area_avail_mu = membrane_info.at['area_avail_mu','Value']
        alpha_numbers_dict = active_membrane_info.alpha_numbers.to_dict()
        enzyme_location = active_membrane_info.location.to_dict()

        membrane_sector = MembraneSector(area_avail_0=area_avail_0,
                                         area_avail_mu=area_avail_mu,
                                         alpha_numbers_dict=alpha_numbers_dict,
                                         enzyme_location=enzyme_location,
                                         max_area=max_membrane_area)

    else:
        membrane_sector = None


    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pamodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_info,
                      translational_sector=translation_enzyme_info,
                       unused_sector=unused_enzyme_info,
                      membrane_sector=membrane_sector,
                      sensitivity=sensitivity, configuration = config
                      )
    return pamodel

def increase_kcats_in_parameter_file(kcat_increase_factor: int,
                                     pam_info_file_path_ori: str,
                                     pam_info_file_path_out: str = None):

    if pam_info_file_path_out is None:
        pam_info_file_path_out= f'{pam_info_file_path_ori.split(".")[0]}_multi.xlsx'

    pam_info_dfs = pd.read_excel(pam_info_file_path_ori, sheet_name=None)

    pam_info_file_ori = pam_info_dfs.pop('ActiveEnzymes')
    pam_info_file_new = pam_info_file_ori.rename({'kcat_values': 'kcat_values_ori'})
    pam_info_file_new['kcat_values'] = pam_info_file_ori['kcat_values']*kcat_increase_factor

    write_mode = 'w'
    kwargs = {}
    if os.path.isfile(pam_info_file_path_out):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(pam_info_file_path_out, engine = 'openpyxl',
                        mode = write_mode, **kwargs) as writer:
        pam_info_file_new.to_excel(writer, sheet_name = 'ActiveEnzymes', index=False)
        for sheet_name, df in pam_info_dfs.items():
            df.to_excel(writer, sheet_name = sheet_name, index=False)