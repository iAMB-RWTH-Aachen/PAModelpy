import pandas as pd
import numpy as np
import cobra
from typing import TypedDict, Literal, Union, Tuple, Iterable, List, Dict, Optional, Type, Any
import re
import os
import ast
from pathlib import Path

from collections import defaultdict
from dataclasses import dataclass, field


from ..PAModel import PAModel
from ..EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from ..configuration import Config

DEFAULT_MOLMASS = 39959.4825 #kDa
DEFAULT_KCAT = 13.7 #s-1, from Bar-Evan et al (2011), the median from BRENDA

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

# def parse_gpr_information(gpr_info:str,
#                           genes:list=None,
#                           enzyme_id:str=None,
#                           gene2protein: dict[str, str] = None,
#                           convert_to_complexes= False) -> tuple[list,list]:
#     #filter out nan entries
#     if not isinstance(gpr_info, str) or isinstance(enzyme_id, float):
#         return None, None
#
#     # #only get the genes associated with this enzyme
#     gpr_list = _parse_gpr(gpr_info)
#     if genes is None: return gpr_list
#
#     gpr_list = _filter_sublists(gpr_list, genes)
#
#     #convert the genes to the associated proteins for a database with single protein names
#     enzyme_relations = []
#     if convert_to_complexes:
#         for sublist in gpr_list:
#             enz_sublist = []
#             for item in sublist:
#                 if item in gene2protein.keys():
#                     if '_' not in gene2protein[item]:
#                         enz_sublist.append(gene2protein[item])
#                     elif gene2protein[item].split('_') not in enzyme_relations:
#                         enzyme_relations += gene2protein[item].split('_')
#
#             enzyme_relations += [enz_sublist]
#     # convert genes to proteins for building a model
#     elif any([len(info)>1 for info in gpr_list]):
#         enzyme_relations = [enzyme_id.split('_')]
#     else:
#         enzyme_relations = [[enzyme_id]]
#         # enzyme_relations = _filter_sublists(enzyme_relations, enzyme_id.split('_'), how='all')
#     return sorted(gpr_list), sorted(enzyme_relations)

def parse_gpr_information(
    gpr_info: str,
    genes: Optional[List[str]] = None,
    enzyme_id: Optional[str] = None,
    gene2protein: Optional[Dict[str, str]] = None,
    convert_to_complexes: bool = False
) -> Tuple[List[List[str]], List[List[str]]]:
    """
    Parses Gene-Protein-Reaction (GPR) information and maps genes to their associated proteins or enzyme complexes.

    Args:
    gpr_info : str
        The GPR string representing the gene associations.
    genes : list[str], optional
        A list of genes to filter the GPR associations.
    enzyme_id : str, optional
        The enzyme ID associated with this reaction.
    gene2protein : dict[str, str], optional
        A dictionary mapping genes to proteins.
    convert_to_complexes : bool, default=False
        Whether to convert genes to protein complexes.

    Returns:
    tuple[list[list[str]], list[list[str]]]
        - A sorted list of gene associations (GPR structure).
        - A sorted list of enzyme relations (protein mappings).
    """

    # Filter out invalid entries
    if not isinstance(gpr_info, str) or isinstance(enzyme_id, float):
        return None, None

    # Step 1: Parse the GPR string
    gpr_list = _parse_gpr(gpr_info)

    # Step 2: Filter genes if provided
    if genes is not None:
        gpr_list = _filter_sublists(gpr_list, genes)

    # Step 3: Convert genes to associated proteins
    enzyme_relations = _map_genes_to_proteins(gpr_list, gene2protein, enzyme_id, convert_to_complexes)

    return sorted(gpr_list), sorted(enzyme_relations)


def _map_genes_to_proteins(
    gpr_list: List[List[str]],
    gene2protein: Optional[Dict[str, str]],
    enzyme_id: Optional[str],
    convert_to_complexes: bool
) -> List[List[str]]:
    """
    Converts genes in a GPR list to their associated proteins or enzyme complexes.

    Args:
    gpr_list : list[list[str]]
        Parsed GPR information as a nested list.
    gene2protein : dict[str, str], optional
        Dictionary mapping genes to proteins.
    enzyme_id : str, optional
        The enzyme ID for this reaction.
    convert_to_complexes : bool
        Whether to convert genes to enzyme complexes.

    Returns:
    list[list[str]]
        A nested list representing the enzyme relationships.
    """

    if gene2protein is None:
        gene2protein = {}

    if convert_to_complexes:
        enzyme_relations = [
            [gene2protein[item] for item in sublist if item in gene2protein]
            for sublist in gpr_list
        ]

    elif any(len(info) > 1 for info in gpr_list):  # Complex enzymes
        enzyme_relations = [enzyme_id.split("_")]
    else:
        enzyme_relations = [[enzyme_id]]

    return enzyme_relations


def get_protein_gene_mapping(enzyme_db: pd.DataFrame,
                             model) -> tuple[dict, dict]:
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
            if isinstance(gene_id, str): #there is a gene id associated
                gene2protein[gene_id] = enzyme_id if pd.notna(enzyme_id) else f"Enzyme_{gene_id}"
            else: # if there is no gene associated, use the model genes
                for gene in rxn.genes:
                    enz_id = enzyme_id if pd.notna(enzyme_id) else f"Enzyme_{gene.id}"
                    gene2protein[gene.id] = enz_id
                    protein2gene[enz_id].append(gene.id)

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

def _filter_sublists(nested_list:list[list[str]], target_list:list[str],
                     how: Literal['any', 'all'] = 'any') -> list[str]:
    """
    Keeps only sublists from `nested_list` that contain at least one identifier present in `target_list`.

    Args:
        nested_list (list of lists): The main list of sublists to filter.
        target_list (list of str): The flat list of identifiers to include.

    Returns:
        list of lists: The filtered nested list.
    """
    if how == 'any':
        filtered_list = [sublist for sublist in nested_list if any([item in target_list for item in sublist])]
    else:
        filtered_list = [sublist for sublist in nested_list if all([item in target_list for item in sublist])]

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
                                                          protein2gpr:defaultdict[str,list],
                                                          gene2protein:dict[str, str]
                                                          ) -> Tuple[dict[str, ReactionInformation],defaultdict[str,list]]:
    for rxn in model.reactions:
        rxn_id = _extract_reaction_id(
            rxn.id)  # some reactions ids are associated with copy numbers, only filter for the actual reaction id
        rxn_genes = [g.id for g in rxn.genes if 's0001' not in g.id]
        if (
                rxn_id in rxn_info2protein.keys()
                or 'EX'.lower() in rxn.id.lower()  # is the reaction an exchange with the environment?
                or 'BIOMASS' in rxn.id  # is the reaction a pseudoreaction?
                or len(rxn_genes) == 0 # is the reaction associated with enzymes?
                or rxn_id == 'ATPM' # is the reaction associated to the ATP maintenance pseudoreaction?
        ): continue

        print('No enzyme information found for reaction: ' + rxn.id)

        rxn_info = ReactionInformation(rxn.id)
        rxn_info.model_reactions = [rxn]
        if 'tpp' in rxn_id:
            kwargs = {'kcat_f':1e6,'kcat_b':1e6} #make sure transporters do not burden the model
        elif rxn_info.check_reaction_reversibility() > 0:
            kwargs = {'kcat_b':1e6}
        elif rxn_info.check_reaction_reversibility() < 0:
            kwargs = {'kcat_f': 1e6}
        else:
            kwargs = {}


        enzyme_info = enzyme_information(rxn.id, rxn._genes, **kwargs)
        rxn_info.enzymes[enzyme_info['enzyme_id']] = enzyme_info

        gpr_info,_ = parse_gpr_information(str(rxn.gpr),
                                         genes=[g.id for g in rxn._genes],
                                         enzyme_id=enzyme_info['enzyme_id'],
                                         gene2protein=gene2protein)

        if gpr_info is None: gpr_info = [[f'gene_{rxn.id}']]

        rxn_info2protein[rxn.id] = rxn_info
        protein2gpr[enzyme_info['enzyme_id']] = gpr_info

    return rxn_info2protein, protein2gpr

def _get_genes_for_proteins(enzyme_db: pd.DataFrame, model) -> dict:
    protein2gene = defaultdict(list)
    gene2protein = defaultdict(list)
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
            gene2protein[gene].append(enzyme_id)
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

def _order_enzyme_complex_id(enz_id:str,
                              other_enzyme_id_pattern: str = r'E[0-9][0-9]*')->str:
    # Define the regex pattern for protein IDs, obtained from UniProtKB, 2024-08-07
    # https://www.uniprot.org/help/accession_numbers
    protein_id_pattern = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'

    # Define the regex pattern to match protein IDs
    protein_id_regex = re.compile(protein_id_pattern + r'|' + other_enzyme_id_pattern)
    proteins = [obj.group(0).strip('_') for obj in re.finditer(protein_id_regex, enz_id)]
    return "_".join(sorted(proteins))


def parse_reaction2protein(enzyme_db: pd.DataFrame,
                           model: cobra.Model,
                           other_enzyme_id_pattern: str = r'(E\d+[a-z]?|Enzyme_\w+)') -> dict:
    rxn_info2protein = {}
    protein2gpr = defaultdict(list)
    #remove copy number substrings from the reaction to make it matchable to enzyme information
    filtered_model_reactions = [_extract_reaction_id(r.id) for r in model.reactions]

    # replace NaN values with unique identifiers
    enzyme_db.loc[enzyme_db['enzyme_id'].isnull(), 'enzyme_id'] = [f'E{i}' for i in
                                                                   range(enzyme_db['enzyme_id'].isnull().sum())]

    enzyme_db.loc[enzyme_db['gene'].isnull(), 'gene'] = [[f'gene_{i}'] for i in
                                                                   range(enzyme_db['gene'].isnull().sum())]


    #make sure all enzyme complexes have an id ordered in a structured way
    enzyme_db['enzyme_id'] = enzyme_db['enzyme_id'].apply(_order_enzyme_complex_id,
                                                        other_enzyme_id_pattern = other_enzyme_id_pattern)
    protein2gene, gene2protein = _get_genes_for_proteins(enzyme_db, model)

    # parse the information for all gene-protein-reaction relations in the dataframe
    for (rxn_id, enzyme_id, gpr), catalytic_reaction_info in enzyme_db.groupby(
            ['rxn_id', 'enzyme_id', 'GPR']):
        # only parse those reactions which are in the model
        if rxn_id not in filtered_model_reactions: continue

        rxn_info = rxn_info2protein.setdefault(rxn_id, ReactionInformation(rxn_id))
        #sometimes, multiple copies are associated with a single reaction
        rxns = rxn_info.get_reaction_from_model(model)
        #the genes are generally stored in a list, which needs to be recovered from the string formatted column
        genes = catalytic_reaction_info.gene.iloc[0]
        if isinstance(genes, str) and genes[-1] == ']':
            genes = ast.literal_eval(catalytic_reaction_info.gene.iloc[0])
        elif isinstance(genes, str):
            genes = [genes]

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
    rxn_info2protein, protein2gpr = _check_if_all_model_reactions_are_in_rxn_info2protein(model,
                                                                                          rxn_info2protein,
                                                                                          protein2gpr,
                                                                                          gene2protein)

    #convert the dataobject dict to a normal dict for correct parsing by PAModelpy
    rxn2protein = {rxn_id: dict(rxn_info.enzymes) for rxn_id, rxn_info in rxn_info2protein.items()}
    return rxn2protein, dict(protein2gpr)


# Function to parse GPR and determine multimer
def merge_enzyme_complexes(df, gene2protein):
    collapsed_rows = []

    for (rxn_id, direction), group in df.groupby(['rxn_id', 'direction']):
        for _, row in group.iterrows():
            #skip nan entries
            if isinstance(row.GPR, float) or isinstance(row.gene, float):
                continue
            # make sure the enzymes are not counted double by dropping potential duplicated
            group = group.drop_duplicates(['enzyme_id', 'rxn_id', 'direction'])
            #parse GPR
            gpr_list, enzyme_relations = parse_gpr_information(
                row['GPR'], row['gene'], row['enzyme_id'], gene2protein, convert_to_complexes=True
            )
            # Collapse "and" relationships into multimer ID
            if enzyme_relations and not all(all([isinstance(e, float) for e in er]) for er in enzyme_relations):
                for gene_list, enzyme_list in zip(gpr_list, enzyme_relations):
                    row_copy = row.copy()
                    #Replace gene with multimer, make sure it is sorted to avoid duplicate entries of the same enzyme complex
                    row_copy['enzyme_id'] = "_".join(sorted(enzyme_list))
                    # add all the annotations to the corresponding gene
                    row_copy['gene'] = gene_list

                    if len(enzyme_list) > 1:
                        # Compute the sum of molMass and length only if it's a complex (more than one enzyme)
                        molmass_sum = group.loc[
                            df['enzyme_id'].isin(enzyme_list), 'molMass'].sum()
                        length_sum = group.loc[
                            df['enzyme_id'].isin(enzyme_list), 'Length'].sum()
                        kcat_mean = group.loc[
                            df['enzyme_id'].isin(enzyme_list), 'kcat_values'].mean()
                    else:
                        # Keep the original molMass/Length if it's a single enzyme
                        molmass_sum = row['molMass']
                        length_sum = row['Length']
                        kcat_mean = row['kcat_values']
                    row_copy['molMass'] = molmass_sum
                    row_copy['Length'] = length_sum
                    row_copy['kcat_values'] = kcat_mean

                    collapsed_rows.append(row_copy)
            else:
                collapsed_rows.append(row)

    collapsed_df = pd.DataFrame(collapsed_rows).drop_duplicates(['rxn_id', 'enzyme_id', 'direction'])
    return collapsed_df

def build_coarse_grained_sector_object(
    pam_info_file:Union[Path,str],
    sheet_name: str,
    sector_cls: Type,
    config: Any,
    prefix: str,
) -> Optional[Any]:
    """
    Build a coarse-grained protein sector from Excel with at least the following values in the 'Parameter' column:

        <prefix>_mu , <prefix>_0 , mol_mass , id_list

    Args:
        pam_info_file
            Path to the PAM‑info Excel workbook.
        sheet_name
            Excel sheet that holds the parameters.
        sector_cls
            Class to instantiate (e.g. ``TransEnzymeSector``).
        config
            The configuration object to set model-specific parameters
        prefix
            prefix of parameters required to build the sector (``'tps'`` for Translational,
            ``'ups'`` for UnusedEnzyme, …)

    Returns:
        An instance of ``sector_cls`` or ``None`` if the sheet could not be read
        or a required key is missing.
    """
    try:
        df = pd.read_excel(pam_info_file, sheet_name=sheet_name).set_index("Parameter")
    except Exception as exc:
        raise KeyError(f"Could not read sheet '{sheet_name}': {exc}")

    slope_key   = f"{prefix}_mu"
    intercept_key = f"{prefix}_0"

    # Helper that works for both index‑based and raw‑DataFrame access
    def _get(key: str):
        return df.at[key, "Value"]

    try:
        mu_val   = _get(slope_key)
        zero_val = _get(intercept_key)
        mol_val  = _get("mol_mass")
        id_val   = _get("id_list")
    except Exception as exc:
        raise KeyError(f"Missing required key while reading '{sheet_name}': {exc}")

    kw = {
        "id_list": [id_val],
        f"{prefix}_mu": [mu_val],
        f"{prefix}_0": [zero_val],
        "mol_mass": [mol_val],
        "configuration": config,
    }

    try:
        return sector_cls(**kw)
    except Exception as exc:
        raise Exception(f"Could not instantiate {sector_cls.__name__}: {exc}")

def build_active_enzyme_sector_from_pam_info_file(pam_info_file:Union[Path,str],
                                                  model: 'Model',
                                                  config: Config,
                                                  adjust_reaction_ids:Optional[bool] = False,
                                                  enzyme_db:Optional[pd.DataFrame]=None,) -> ActiveEnzymeSector:
    if enzyme_db is None:
        enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')
        # for some models, the reaction ids should not include 'pp' or 'ex'
        if adjust_reaction_ids:
            enzyme_db['rxn_id'] = enzyme_db['rxn_id'].apply(_check_rxn_identifier_format)
    # create enzyme objects for each gene-associated reaction
    rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, model)

    return ActiveEnzymeSector(rxn2protein=rxn2protein, protein2gene=protein2gene,
                                            configuration=config)

def set_up_pam(pam_info_file:Union[Path,str] = '',
               model:Union[str, cobra.Model] = 'Models/iML1515.xml',
               config:Config = None,
               total_protein: Union[bool, float] = True,
               active_enzymes: bool = True,
               translational_enzymes: bool = True,
               unused_enzymes: bool = True,
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

    active_enzyme_info = build_active_enzyme_sector_from_pam_info_file(pam_info_file = pam_info_file,
                                                                           model = model,
                                                                           adjust_reaction_ids = adjust_reaction_ids,
                                                                           enzyme_db = enzyme_db,
                                                                           config = config) if active_enzymes else None

    translation_enzyme_info = build_coarse_grained_sector_object(pam_info_file = pam_info_file,
                                           sheet_name='Translational',
                                           sector_cls=TransEnzymeSector,
                                           prefix = 'tps',
                                           config=config) if translational_enzymes else None


    unused_enzyme_info = build_coarse_grained_sector_object(pam_info_file = pam_info_file,
                                           sheet_name='UnusedEnzyme',
                                           sector_cls=UnusedEnzymeSector,
                                           prefix = 'ups',
                                           config=config) if unused_enzymes else None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pamodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_info,
                      translational_sector=translation_enzyme_info,
                       unused_sector=unused_enzyme_info,
                      sensitivity=sensitivity, configuration = config
                      )
    return pamodel

def set_up_core_pam(pam_info_file:str = '',
               model:Union[str, cobra.Model] = 'Models/e_coli_core.json',
               config:Config = None,
               total_protein: Union[bool, float] = True,
               active_enzymes: bool = True,
               translational_enzymes: bool = True,
               unused_enzymes: bool = True,
               sensitivity:bool = True,
               enzyme_db:pd.DataFrame = None,
               adjust_reaction_ids:bool = True) -> PAModel:


    if config is None:
        config = Config()
        config.reset()

    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    config.BIOMASS_REACTION = BIOMASS_REACTION
    TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]

    #setup model if a path is provided
    if isinstance(model, str):
        model = cobra.io.load_json_model(model)

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
        # translational protein sector parameter (substrate dependent)
        id_list_tps = ['EX_glc__D_e']
        tps_0 = [0.04992]  # g/gDW
        tps_mu = [-0.002944]  # g h/gDW -> transformed to match glucose uptake variable
        molmass_tps = [405903.94]  # g/mol

        # translational protein sector
        translation_enzyme_info = TransEnzymeSector(
            id_list=id_list_tps,
            tps_0=tps_0,
            tps_mu=tps_mu,
            mol_mass=molmass_tps,
            configuration=config
        )
    else:
        translation_enzyme_info = None

    if unused_enzymes:
        id_list_ups = [BIOMASS_REACTION]
        ups_0 = [0.0407]  # g/gDW
        ups_mu = [-0.0214]  # g h/gDW -> negative relation with growth rate
        molmass_ups = [405903.94]  # g/mol

        unused_enzyme_info = UnusedEnzymeSector(
            id_list=id_list_ups,
            ups_0=ups_0,
            ups_mu=ups_mu,
            mol_mass=molmass_ups,
            configuration=config
        )
    else:
        unused_enzyme_info = None
    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    coremodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_info,
                      translational_sector=translation_enzyme_info,
                       unused_sector=unused_enzyme_info,
                      # membrane_sector=membrane_sector,
                      sensitivity=sensitivity, configuration = config
                      )
    return coremodel

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