from cobra import Object, DictList, Gene, Model
from optlang import Constraint, Variable
from optlang.symbolics import Zero
from typing import Optional, Union
from src.PAModelpy.Enzyme import Enzyme

class Transcript(Object):
    def __init__(self, id: str, gene:Union[Gene, list], enzymes: DictList, length: int):
        DEFAULT_TRANSCRIPT_LENGTH = 750 #nt
        super().__init__(name = id)
        if isinstance(gene, Gene): gene = [gene]
        self.id = id
        self.genes = gene
        self.enzymes = enzymes
        if length>0: self.length = length
        else: self.length = DEFAULT_TRANSCRIPT_LENGTH
        self.mrna_variable = Variable(name = id, lb =0)
        self.f_max = 1
        self.f_min = 0.5
        self._model = None
        self._constraints = {}
        self._lumped_transcripts = DictList()

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model
        self._model.add_cons_vars(self.mrna_variable)
        self._model.transcripts.append(self)
        # check if genes are in the model
        for gene in self.genes:
            if not gene in model.genes: model.genes += self.genes
        # check if enzymes are in the model
        for enzyme in self.enzymes:
            if not enzyme in model.enzymes: model.add_enzyme(enzyme)
            if not isinstance(self, LumpedTranscript):
                transcript_associated_with_enzyme = model.get_transcripts_associated_with_enzyme(enzyme)
                for relation, transcripts in transcript_associated_with_enzyme.items():
                    model.make_mrna_min_max_constraint(enzyme, transcripts, relation)
            else:
                model.make_mrna_min_max_constraint(enzyme, [self])
            # print(self.mrna_variable.lb)
            # for key, value in self._constraints.items():
            #     print(value)

    def change_concentration(self, concentration:float, error: float = 0) -> None:
        """ Setting the concentration of an mRNA species

        Changes the mRNA concentration and thereby influencing the enzyme concentration and reaction rate. The bounds
        on the mrna variable can be made more flexible by including a (measurement) error

        :param concentration:
        :return:
        """
        if len(self._lumped_transcripts)>0:
            # if there is an and relationship, the lumped transcript should be adjusted instead of the individual transcript
            for lumped_transcript in self._lumped_transcripts:
                if lumped_transcript.mrna_variable.lb > concentration-error: #only if this transcript is the limiting transcript
                    print(f'changing lumped transcript concentration from {lumped_transcript.mrna_variable.lb} to {concentration-error}')
                    lumped_transcript.mrna_variable.ub = concentration + error
                    lumped_transcript.mrna_variable.lb = concentration - error
        else:
            self.mrna_variable.ub = concentration+error
            self.mrna_variable.lb = concentration-error
        for enzyme in self.enzymes:
            self._model.correct_unused_enzymes_for_measured_mrna(enzyme)
        self._model.solver.update()

    def reset_concentration(self):
        if len(self._lumped_transcripts)>0:
            # if there is an and relationship, the lumped transcript should be adjusted instead of the individual transcript
            for lumped_transcript in self._lumped_transcripts:
                lumped_transcript.mrna_variable.ub = 1e6
                lumped_transcript.mrna_variable.lb = -1e6
        else:
            self.mrna_variable.ub = 1e6
            self.mrna_variable.lb = -1e6
        self._model.solver.update()

class LumpedTranscript(Transcript):
    def __init__(self, id: str, gene:Union[Gene, list], enzymes: DictList, length: int):

        super().__init__(id = id,
                         gene = gene,
                         enzymes =enzymes,
                         length = length)
        self._transcripts = DictList()