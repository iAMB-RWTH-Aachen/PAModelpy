import pytest

from src.PAModelpy import CatalyticEvent

def test_if_catalytic_events_splits_rxns_from_enzymes():
    # Arrange
    sut = CatalyticEvent
    catalytic_reaction_ids = [
    "CYTBO3_4pp", "CE_CYTBO3_4pp_P0ABJ1_P0ABJ5_P0ABJ7_P0ABJ8", "CE_CYTBO3_P18200_P0AEZ7",
    "CE_ME_P23N45", "ME",
    "CE_REACTID_O123A4", "CE_REACTID_A1A2A3",
    "LPLIPAL2ATE140","CE_LPLIPAL2ATE140_P23N45"
    ]
    rxn_ids = ["CYTBO3_4pp", "CYTBO3_4pp", "CYTBO3",
               "ME", "ME",
               "REACTID", "REACTID",
               "LPLIPAL2ATE140", "LPLIPAL2ATE140"]

    # Apply
    parsed_rxn_ids = []
    for cr_id in catalytic_reaction_ids:
        parsed_rxn_ids.append(sut._extract_reaction_id_from_catalytic_reaction_id(cr_id))

    # Assert
    assert all([rxn_id_sut == rxn_id for rxn_id_sut, rxn_id in zip(parsed_rxn_ids, rxn_ids)])