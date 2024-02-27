from Scripts.tam_generation import set_up_toy_tam, set_up_ecolicore_tam
from Scripts.pam_generation import set_up_ecolicore_pam, set_up_toy_pam


def test_set_up_toy_tam_works():
    sut = set_up_toy_tam

def test_set_up_ecolicore_tam_works():
    sut = set_up_ecolicore_tam()
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

    # Arrange
    assert 'optimal' == toy_pam.solver.status
    assert 'optimal' == sut.solver.status
    assert toy_pam.objective.value == sut.objective.value

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
    for enzyme in ecolicore_pam.enzyme_variables:
        print(enzyme.id, enzyme.flux)

    # Arrange
    assert 'optimal' == ecolicore_pam.solver.status
    assert 'optimal' == sut.solver.status
    assert ecolicore_pam.objective.value == sut.objective.value
