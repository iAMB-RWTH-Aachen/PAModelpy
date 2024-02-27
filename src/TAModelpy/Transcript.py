from cobra import Object, DictList, Gene, Model
from optlang import Constraint, Variable
from optlang.symbolics import Zero
from typing import Optional
from src.PAModelpy.Enzyme import Enzyme

class Transcript(Object):
    def __init__(self, id: str, gene:Gene, enzymes: DictList, length: int):

        super().__init__(name = id)
        self.id = id
        self.gene = gene
        self.enzymes = enzymes
        self.length = length
        self.mrna_variable = Variable(name = id, lb =0)
        self.f_max = 1
        self.f_min = 0.5
        self._model = None
        self._constraints = {}

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model
        self._model.add_cons_vars(self.mrna_variable)
        # check if gene is in the model
        if not self.gene in model.genes: model.genes.append(self.gene)
        # check if enzymes are in the model
        for enzyme in self.enzymes:
            if not enzyme in model.enzymes: model.add_enzyme(enzyme)
            # #check the gpr relations: is an additional gene required or is there 'competition' with another gene?
            # other_mrnas= model.get_transcripts_associated_with_enzyme(enzyme)
            # #TODO how to handle other mrna relations?
            # # add the mRNA relations to the model
            # self._model.make_mrna_min_max_constraint(enzyme, self)

    def change_concentration(self, concentration:float, error: float = 0) -> None:
        """ Setting the concentration of an mRNA species

        Changes the mRNA concentration and thereby influencing the enzyme concentration and reaction rate. The bounds
        on the mrna variable can be made more flexible by including a (measurement) error

        :param concentration:
        :return:
        """
        self.mrna_variable.ub = concentration+error
        self.mrna_variable.lb = concentration-error
        self._model.solver.update()
