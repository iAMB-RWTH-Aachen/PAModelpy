from cobra import Reaction, Metabolite
from typing import Optional


class Constraint(Metabolite):
    """
    Class for information about a Constraint in a protein Sector.

    Constraint is a class for holding information similar to
    a metabolite in a cobra.Reaction object.

    Parameters:
        id (str): The identifier to associate with the constraint.
        name (str): A human-readable name.
    """

    # noinspection PyShadowingBuiltins
    def __init__(
        self,
        id: Optional[str] = None,
        formula: Optional[str] = None,
        name: Optional[str] = "",
        charge: Optional[float] = None,
        compartment: Optional[str] = "cytosol",
    ) -> None:
        super().__init__(
            id=id,
            formula=formula,
            name=name,
            charge=charge,
            compartment=compartment,
        )
        self.annotation = {"type": "Constraint"}
