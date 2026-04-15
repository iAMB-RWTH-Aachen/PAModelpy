import os
import cobra.io
import numpy as np
import pandas as pd
import pytest
from typing import List, Dict

from Scripts.toy_ec_pam import build_toy_gem
from src.PAModelpy import PAModel, TransEnzymeSector, UnusedEnzymeSector
from src.PAModelpy.utils.pam_generation import parse_reaction2protein, set_up_pam, merge_enzyme_complexes, build_coarse_grained_sector_object

@pytest.fixture
def translational_df():
    return pd.DataFrame(
        {
            "Parameter": ["id_list", "tps_mu", "tps_0", "mol_mass"],
            "Value": ["T001", 0.1, 5.0, 55000],
        }
    )

@pytest.fixture
def unused_df():
    return pd.DataFrame(
        {
            "Parameter": ["id_list", "ups_mu", "ups_0", "mol_mass"],
            "Value": ["U001", 0.05, 2.0, 50000],
        }
    )
@pytest.fixture
def pam_info_file(tmp_path, translational_df, unused_df) -> str:
    path = tmp_path / "pam_info.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        translational_df.to_excel(writer, sheet_name="Translational", index=False)
        unused_df.to_excel(writer, sheet_name="UnusedEnzyme", index=False)
    return str(path)

def test_if_rxn2protein_info_is_correctly_parsed():
    # Arrange
    toy_enzyme_db = pd.DataFrame({'rxn_id':['R1', 'R2', 'R2', 'R3', 'R2'],
                                  'enzyme_id': ['E1', 'E2a', 'E2b_E2c', 'E3', 'E2a'],
                                  'gene': [['gene1'],['gene2a'], ['gene2b', 'gene2c'], ['gene3'], ['gene2a']],
                                  'GPR': ['gene1','gene2a or (gene2b and gene2c)', 'gene2a or (gene2b and gene2c)',
                                          'gene3', 'gene2a or (gene2b and gene2c)'],
                                  'molMass':[1, np.nan,2,3,2],
                                  'kcat_values':[1,2,np.nan,3,1.5],
                                  'direction':['f','f', 'f', 'f', 'b']
                                  }
                                 )
    enzyme_id_pattern = r'E\d+[a-z]?'
    toy_model = build_toy_gem()

    expected_rxn2protein = {
        'R1':  {
            'E1': {
                'enzyme_id': 'E1', 'f': 1.0, 'b': 0,
                'genes': ['gene1'], 'protein_reaction_association': [['E1']],
                'molmass': 1.0}
        },
        'R2': {
            'E2a': {
                'enzyme_id': 'E2a', 'f': 2.0, 'b': 1.5,
                'genes': ['gene2a'], 'protein_reaction_association': [['E2a']],
                'molmass': 39959.4825},
            'E2b_E2c': {
                'enzyme_id': 'E2b_E2c', 'f': 0, 'b': 0,
                'genes': ['gene2b', 'gene2c'], 'protein_reaction_association': [['E2b', 'E2c']],
                'molmass': 2.0}
        },
        'R3': {
            'E3': {
                'enzyme_id': 'E3', 'f': 3.0, 'b': 0,
                'genes': ['gene3'], 'protein_reaction_association': [['E3']],
                'molmass': 3.0}
        }
    }
    expected_protein2gpr = {'E1': [['gene1']], 'E2a': [['gene2a']], 'E2b_E2c': [['gene2b', 'gene2c']], 'E3': [['gene3']]}

    # Apply
    rxn2protein, protein2gpr = parse_reaction2protein(toy_enzyme_db, toy_model,
                                                      other_enzyme_id_pattern = enzyme_id_pattern)
    # Assert
    for output_dict, expected_dict in zip([rxn2protein, protein2gpr], [expected_rxn2protein, expected_protein2gpr]):
        assert all([expected_dict[key] == value for key, value in output_dict.items()])


@pytest.mark.parametrize('path_to_model', [
    (os.path.join('Models', 'e_coli_core.json')),
    (os.path.join('Models', 'iML1515.xml'))
])
def test_if_pam_can_be_build_from_path_to_gem(path_to_model:str):
    pam = PAModel(path_to_model,
                  sensitivity = False)
    pam.change_reaction_bounds('EX_glc__D_e', -10,0)
    pam.optimize()
    assert pam.objective.value > 0


def test_if_set_up_pam_can_build_ecolicore_pam():
    #Arrange
    pam_data_file = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_241209_core.xlsx')
    ecolicore_gem = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))

    #Apply
    ecolicore_pam = set_up_pam(pam_data_file,
                               ecolicore_gem,
                               total_protein = 0.1699,
                               sensitivity=False,
                               adjust_reaction_ids=True)

    ecolicore_pam.optimize()

    #Assert
    assert ecolicore_pam.objective.value > 0

def test_if_set_up_pam_can_build_iML1515():
    #Arrange
    pam_data_file = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_241209.xlsx')
    iml1515 = os.path.join('Models', 'iML1515.xml')

    #Apply
    pam = set_up_pam(pam_data_file,
                               iml1515,
                               sensitivity=False,
                               adjust_reaction_ids=True)

    pam.optimize()

    #Assert
    assert pam.objective.value > 0

def test_if_merge_enzyme_complexes_merges_enzyme_ids():
    # Arrange
    toy_enzyme_db = pd.DataFrame({'rxn_id': ['R1', 'R2', 'R2', 'R2', 'R3'],
                                  'enzyme_id': ['E1', 'E2a', 'E2b', 'E2c', 'E3'],
                                  'gene': [['gene1'], ['gene2a'], ['gene2b', 'gene2c'], ['gene2b', 'gene2c'], ['gene3']],
                                  'GPR': ['gene1', 'gene2a or (gene2b and gene2c)', 'gene2a or (gene2b and gene2c)', 'gene2a or (gene2b and gene2c)',
                                          'gene3'],
                                  'molMass': [1, 2, 2, 2, 3],
                                  'kcat_values': [1, 2, np.nan, 3, 1.5],
                                  'direction': ['f', 'f','f', 'f', 'f'],
                                  'Length': [1]*5
                                  }
                                 )
    gene2protein = {'gene1':'E1','gene2a': 'E2a','gene2b': 'E2b','gene2c': 'E2c', 'gene3':'E3'}

    # Apply
    merged_enzyme_db = merge_enzyme_complexes(toy_enzyme_db, gene2protein)

    # Assert
    assert_enzyme_complexes_are_merged(merged_enzyme_db)
    assert_isozymes_are_not_merged(merged_enzyme_db)
    assert_molmass_for_enzyme_complexes_are_summed(merged_enzyme_db)
    assert_molmass_for_isozymes_are_not(merged_enzyme_db)

@pytest.mark.parametrize('sheetname, gene2protein, enzyme_complexes, len_enzyme_db, molmass_to_check',[
    (
            'mixed_and_or_gpr',
            {'cg3361':'Enzyme_cg3361' ,  'cg3360':'Enzyme_cg3360', 'cg3359':'Enzyme_cg3359'},
            ['Enzyme_cg3359_Enzyme_cg3360', 'Enzyme_cg3359_Enzyme_cg3361', 'Enzyme_cg3361'],
            3,
            {'Enzyme_cg3359':2*39959.4825}
     ),
    (
            'f_b_multiple_complexes_atps',
            {'b3739': 'P0ABC0', 'b3731': 'P0A6E6', 'b3733': 'P0ABA6', 'b3735': 'P0ABA4', 'b3734': 'P0ABB0',
             'b3732': 'P0ABB4', 'b3738': 'P0AB98', 'b3736': 'P0ABA0', 'b3737': 'P68699'},
            ['P0A6E6_P0AB98_P0ABA0_P0ABA4_P0ABA6_P0ABB0_P0ABB4_P0ABC0_P68699','P0A6E6_P0AB98_P0ABA0_P0ABA4_P0ABA6_P0ABB0_P0ABB4_P68699'],
            4,
            {'P0ABC0':240979}
    )
])
def test_if_merge_enzyme_complex_parses_complex_gprs(sheetname:str,
                                                     gene2protein:Dict[str,str],
                                                     enzyme_complexes: List[str],
                                                     len_enzyme_db: int,
                                                     molmass_to_check: Dict[str, float]):
    # Arrange
    toy_enzyme_db = pd.read_excel(os.path.join('tests', 'data', 'enzyme_complexes_to_merge_info.xlsx'),
                                  sheet_name=sheetname
                                  )
    print(toy_enzyme_db[['rxn_id', 'enzyme_id', 'direction', 'molMass', 'Length']].to_markdown())
    # Apply
    merged_enzyme_db = merge_enzyme_complexes(toy_enzyme_db, gene2protein)
    print(merged_enzyme_db[['rxn_id', 'enzyme_id', 'direction', 'molMass', 'Length']].to_markdown())
    # Assert
    for complex in enzyme_complexes:
        assert (complex==merged_enzyme_db.enzyme_id).any()
    assert len(merged_enzyme_db) == len_enzyme_db
    for enzyme_id, molmass in molmass_to_check.items():
        assert (merged_enzyme_db.molMass[merged_enzyme_db.enzyme_id.str.contains(enzyme_id)] == molmass).all()

@pytest.mark.parametrize('sheet_name,sector_object, prefix', [
    ("Translational", TransEnzymeSector,"tps"),
    ("UnusedEnzyme", UnusedEnzymeSector, "ups")
])
def test_build_translational_sector_success(pam_info_file, sheet_name, sector_object, prefix, translational_df, unused_df):
    sector = build_coarse_grained_sector_object(
        pam_info_file=pam_info_file,
        sheet_name=sheet_name,
        sector_cls=sector_object,
        config={"foo": "bar"},
        prefix=prefix,
    )
    assert isinstance(sector, sector_object)

    # –‑> its stored kwargs must match exactly what the function builds
    expected = translational_df if sheet_name=="Translational" else unused_df
    for i, row in expected.iterrows():
        attr = getattr(sector, row.Parameter)
        assert row["Value"] == attr if not isinstance(attr, list) else attr[0]



#########################################################################################################################
# HELPER FUNCTION
#########################################################################################################################
def assert_enzyme_complexes_are_merged(enzyme_db: pd.DataFrame):
    assert ('E2b_E2c'==enzyme_db.enzyme_id).any()
    assert not ('E2b'==enzyme_db.enzyme_id).any()
    assert not ('E2c'==enzyme_db.enzyme_id).any()

def assert_isozymes_are_not_merged(enzyme_db: pd.DataFrame):
    assert ('E2a'==enzyme_db.enzyme_id).any()
    assert not enzyme_db.enzyme_id.str.contains('_E2a', na=False).any()
    assert not enzyme_db.enzyme_id.str.contains('E2a_', na=False).any()

def assert_molmass_for_enzyme_complexes_are_summed(enzyme_db: pd.DataFrame):
    assert (enzyme_db.molMass[enzyme_db.enzyme_id=='E2b_E2c'] == 4).all()

def assert_molmass_for_isozymes_are_not(enzyme_db: pd.DataFrame):
    assert (enzyme_db.molMass[enzyme_db.enzyme_id=='E2a'] == 2).all()


