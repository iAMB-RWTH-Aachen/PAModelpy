from cobra import Object, Gene, DictList

class Transcript(Object):

    def __init__(self,
                 id:str,
                 gene: Gene,
                 enzymes: DictList,
                 length: int):
        self.id = id
        self.gene = gene
        self.enzymes = enzymes
        self.length = length