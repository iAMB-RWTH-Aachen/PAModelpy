"""Define objects inherited from the cobra.Object class"""

from cobra import Object, DictList
from typing import Optional



class EnzymaticReaction(Object):
    """The EnzymaticReaction is used in a cobra.model to
    get enzymatic information for a specific reaction using the reaction ID
    """

    def __init__(
            self,
            id: Optional[str] = None, # reaction id
            name: str = ""
        ) -> None:
        # initialize the EnzymaticReaction
        super().__init__(id, name)

        # initialize the EnzymaticReactionObject attributes
        self.enzymatic_reactions = DictList()
        self.enzymes = DictList()

    @property
    def enzyme_concentration(self) -> float:
        """returns the total of enzymes actively catalyzing the reaction
        Only the part of the enzymes are considered which participate in the reaction
        Returns
        -------
        float
            Enzyme concentration
        """
        concentration = 0.0
        for enzymatic_reaction in self.enzymatic_reactions:
            concentration += enzymatic_reaction.concentration

        return concentration







