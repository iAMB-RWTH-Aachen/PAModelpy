from typing import List, Optional, Union, Dict, Iterable
from cobra import Object, DictList, Gene

from src.TAModelpy.Transcript import Transcript
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector, CustomSector, Sector
from src.PAModelpy.configuration import Config

class TAModel(PAModel):
    def __init__(self, id_or_model: Union[str, "Model", None] = None,
                 name: Optional[str] = None,
                 p_tot: Optional[float] = Config.P_TOT_DEFAULT,
                 sensitivity: bool = True,
                 active_sector: Optional[ActiveEnzymeSector] = None,
                 translational_sector: Optional[TransEnzymeSector] = None,
                 unused_sector: Optional[UnusedEnzymeSector] = None,
                 custom_sectors: Union[List, CustomSector] = None,
                 gene2transcript: dict = None,
                 configuration=Config):

        self.transcripts = DictList

        if isinstance(id_or_model, PAModel):
            self.build_tamodel_from_pamodel(id_or_model, gene2transcript)
        else:
            super().__init__(id_or_model= id_or_model,
                     name= name,
                     p_tot= p_tot,
                     sensitivity =sensitivity,
                     active_sector = active_sector,
                     translational_sector = translational_sector,
                     unused_sector = unused_sector,
                     custom_sectors=custom_sectors,
                     configuration = configuration)

            self.add_transcript_information(gene2transcript)

#need to connect enzymes to genes

    def build_tamodel_from_pamodel(self, pamodel: PAModel, gene2transcript: dict) -> None:
        #TODO how to initialize from PAM?
        self.add_transcript_information(gene2transcript)

    def add_transcript_information(self, gene2transcript: dict) -> None:
        for gene_id, transcript_info in gene2transcript.items():
            if not gene_id in self.genes: self.genes.append(Gene(gene_id))
            gene_object = self.genes.get_by_id(gene_id)
            enzymes = self.get_enzymes_by_gene_id(gene_id)
            transcript_object = Transcript(id = transcript_info['id'],
                                           gene = gene_object,
                                           enzymes = enzymes,
                                           length = transcript_info['length'])
            self.transcripts.append(transcript_object)

    def get_enzymes_by_gene_id(self, gene_id: str):
        return [enzyme for enzyme in self.enzymes if gene_id in enzyme.genes]


    def add_ta_variable(self):
        pass

    def add_ta_ub_lb(self):
        pass


class Transcript(Object):

    def __init__(self):
        self.gene = None
        self.reactions = DictList()
        self.enzymes = DictList()
        self.enzyme_variables = DictList()
        self.sequence_length = None
