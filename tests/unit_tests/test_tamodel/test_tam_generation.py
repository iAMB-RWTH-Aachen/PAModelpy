import pytest

from Scripts.tam_generation import set_up_toy_tam, set_up_ecolicore_tam, parse_gpr_information_for_protein2genes, set_up_ecoli_tam
from Scripts.pam_generation import set_up_ecolicore_pam, set_up_toy_pam, set_up_ecoli_pam


def test_gpr_information_is_parsed_correctly():
    # Arrange
    enzyme_id = '1.1.5.3'
    gpr_string ='(b0902 and b0903) or (b0902 and b3114) or (b3951 and b3952) or ((b0902 and b0903) and b2579)'

    parsed_gpr_reference = [['b0902','b0903'], ['b0902','b3114'], ['b3951' ,'b3952'],['b0902', 'b0903','b2579']]

    # Act
    gpr_parsed = parse_gpr_information_for_protein2genes(gpr_string, enzyme_id)

    # Assert
    assert gpr_parsed == parsed_gpr_reference

def test_set_up_toy_tam_works():
    sut = set_up_toy_tam()
    assert True

def test_set_up_ecolicore_tam_works():
    sut = set_up_ecolicore_tam()
    assert True

def test_set_up_ecoli_tam_works():
    sut = set_up_ecoli_tam()
    assert True

def test_toy_tam_without_constraints_has_same_result_as_pam():
    # Arrange
    sut = set_up_toy_tam()
    toy_pam = set_up_toy_pam()

    sut.change_reaction_bounds(rxn_id='R1',
                                   lower_bound=0, upper_bound=1e-3)
    toy_pam.change_reaction_bounds(rxn_id='R1',
                                   lower_bound=0, upper_bound=1e-3)
    # Act
    sut.optimize()
    toy_pam.optimize()

    # Assert
    assert_solution_pam_equals_tam(toy_pam, sut)


def test_ecolicore_tam_without_constraints_has_same_result_as_pam():
    # Arrange
    sut = set_up_ecolicore_tam()
    ecolicore_pam = set_up_ecolicore_pam()

    sut.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                   lower_bound=-10, upper_bound=0)
    ecolicore_pam.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                   lower_bound=-10, upper_bound=0)

    # Act
    sut.optimize()
    ecolicore_pam.optimize()
    # Arrange
    assert_solution_pam_equals_tam(ecolicore_pam, sut)

def test_ecoli_tam_without_constraints_has_same_result_as_pam():
    # Arrange
    sut = set_up_ecoli_tam()
    ecoli_pam = set_up_ecoli_pam()

    sut.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                   lower_bound=-10, upper_bound=0)
    ecoli_pam.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                   lower_bound=-10, upper_bound=0)

    # Act
    sut.optimize()
    ecoli_pam.optimize()

    print('in tam but not in pam')
    for enz in sut.enzymes:
        if enz.id not in ecoli_pam.enzymes:
            print(enz.id)

    print('in pam but not in tam')
    for enz in ecoli_pam.enzymes:
        if enz.id not in sut.enzymes:
            print(enz.id)

    # Arrange
    assert_solution_pam_equals_tam(ecoli_pam, sut)


############################################################################################################################
# HELPER FUNCTIONS
############################################################################################################################

def assert_solution_pam_equals_tam(pam, tam):
    assert 'optimal' == pam.solver.status
    assert 'optimal' == tam.solver.status
    assert pam.objective.value == pytest.approx(tam.objective.value, abs=1e-2)

