"""Upper level Enzyme class for handling enzyme objects in cobra models
Enzyme promiscuity is represented by multiple enzyme objects (enzymatic reactions) per enzyme
"""
import cobra.core
import sys

import optlang
from optlang.symbolics import Zero

sys.path.append('../Package/')
from Variables import CatalyticEvent, EnzymeVariable
from typing import Dict, Union, Optional
from cobra import DictList
from warnings import warn



class Enzyme():
    """Upper level Enzyme object containing information about the enzyme
    and link to the enzyme variables (CatalyticEvents) for each reaction the
    enzyme catalyzes.
    This class is used to generate enzyme instances from kcat values and contains
    information about the forward as well as the backward catalysis.

    There are two scenarios:
    - Promiscuous enzymes: a single enzyme can catalyze multiple reactions
    - Other: a single enzyme catalyzes a single reaction

    Parameters
    -------
    id : str
        Identifier for the enzyme (e.g. Uniprot ID)
    rxn2kcat: Dict
        Dictionary with reaction ID, kcat value pairs for the forward (f)
        and backward (b) reaction
        (Example: {'PGI': {'f': 30, 'b': 0.1}})
    upper_bound: float
        Upper bound for the enzyme variable (default 1000.0)
    name: str
        Name of the enzyme (default None)
    molmass: float
        Molar mass of the enzyme (default 3.947778784340140e04)

    """
    
    
    # constant parameters
    DEFAULT_ENZYME_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]
    
    def __init__(
        self,
        id: str,
        rxn2kcat: Dict,
        upper_bound: Union[int, float] = 1000.0,
        name: Optional[str] = None,
        molmass: Union[int, float] = DEFAULT_ENZYME_MOL_MASS,
    ):

        self.rxn2kcat = rxn2kcat
        self.molmass = molmass
        self.id = id # use e.g. Uniprot ID
        self.name = name
        self.upper_bound = upper_bound
        self.catalytic_event_id = 'CE_' + '{0}' # generic template for the catalytic events IDs
        # initialize CatalyticEvent variables for each associated reaction (promiscuity)
        self.catalytic_events = DictList()
        self.enzyme_variable = None
        self.create_enzyme_variable()
        for rxn_id, kcats in rxn2kcat.items():
            self.create_catalytic_event(rxn_id, kcats)
        self._constraints = {} # dict with constraint_id:optlang.Constraint, key:value pairs.
        self._model = None
        self.enzyme_complex = [] #is the enzyme in a complex?
        self.annotation = {'type':'Constraint'}#you can add an annotation for an enzyme
     
    @property
    def kcat_values(self):
        """returns a dictionary with kcat values for each associated reaction
        """
        return self.get_kcat_values()

    @property
    def concentration(self, units:str = 'mmol/gDW', return_units:bool = False) -> float:
        """returns the enzyme's total concentration
        Any associated reaction is considered

        Parameters
        -------
        units: optional, string
            units in which the concentration is calculated (default is mmol/gDW), other option is 'g/gDW'
        return_units: optional, bool
            determines wheter the units should be returned as well

        Returns
        -------
        float
            Enzyme concentration
        """

        # sum up concentrations (aka fluxes) of all enzyme objects
        concentration = 0.0
        for catalytic_event in self.catalytic_events:
            concentration += catalytic_event.flux()
        if units == 'g/gDW':
            #converting mmol to grams of protein:
            # [g] = [mmol]* 1e-3 [mol/mmol] * MW[g/mol]
            concentration = concentration * 1e-3 * self.molmass
        if return_units:
            return concentration, units
        return concentration

    def create_constraint(self, extension: str = None):
        if extension is None: extension = ''
        else: extension = '_' + extension
        # return cobra.core.Metabolite(id = f'EC_{self.id}{extension}', compartment='Enzymes')
        return optlang.Constraint(Zero, name= f'EC_{self.id}{extension}', lb=0, ub=0)

    def add_catalytic_event(self, ce: CatalyticEvent, kcats: Dict):
        """
        Adding catalytic event associated to a reaction to an enzyme
        Parameters
        ----------
        ce: PAModelpy.Variables.CatalyticEvent
            The catalytic event object to which the enzyme should be added
        kcats: dict
            A list with dicts containing direction, kcat key value pairs
        Returns
        -------

        """
        self.catalytic_events += [ce]
        self.enzyme_variable.add_catalytic_events([ce],[kcats])

    def create_catalytic_event(self, rxn_id: str, kcats: Dict):
        """creates enzyme variables that link to reactions

        Parameters
        ----------
        rxn_id : str
            ID of the associated reaction in the model
        kcats : Dict
            kcat values for the forward and backward reaction
        Returns
        -------
        Variables.CatalyticEvent
            Enzyme variable object
        """
        
        # create unique enzyme object name and id
        catalytic_event_id = self.catalytic_event_id.format(rxn_id)
        if self.name is not None:
            enzyme_object_name = rxn_id + '_' + self.name 
        else:
            enzyme_object_name = rxn_id + '_'  + self.id
        
        # create enzymatic reaction object inherited from cobra.Reaction
        catalytic_event= CatalyticEvent(
            id=catalytic_event_id,
            rxn_id=rxn_id,
            kcats2enzymes={self: kcats},
            name=enzyme_object_name
        )

        self.add_catalytic_event(catalytic_event, kcats)

    def create_enzyme_variable(self):
        """creates enzyme variables that link  enzyme to reactions
        """
        # create enzymatic reaction object inherited from cobra.Reaction
        enzyme_variable = EnzymeVariable(
            id=self.id,
            kcats2rxns=self.rxn2kcat,
            upper_bound=self.upper_bound,
        )

        self.enzyme_variable = enzyme_variable

    def change_kcat_values(self, rxn2kcat: Dict):
        """changes the kcat values for the enzyme
        and updates the enzyme variable (enzymatic reaction) accordingly

        Parameters
        ----------
        rxn2kcat : Dict
            Dictionary with reaction ID, kcat value pairs for the forward (f)
            and backward (b) reaction
            (Example: {'PGI': {'f': 30, 'b': 0.1}})
        """

        # update the enzyme variables
        for rxn_id, kcats in rxn2kcat.items():
            catalytic_event_id = self.catalytic_event_id.format(rxn_id)
            # is there already a link between enzyme and reaction?
            if catalytic_event_id not in self.catalytic_events:
                warn(f'Reaction {rxn_id} is not associated with enzyme {self.id}. Skip')
            else:
                # change kcat values of existing enzyme variable
                self.catalytic_events.get_by_id(catalytic_event_id).change_kcat_values({self: kcats})


    def get_kcat_values(self, rxn_ids: Union[str, list] = None) -> Dict:
        """returns the kcat values for a specific enzyme and all
        enzyme-associated reactions

        Parameters
        ----------
        rxn_ids : str or list
            ID of the reactions for which the kcat values should be returned

        Returns
        -------
        Dict
            kcat values for the forward and backward reaction
        """
        
        if isinstance(rxn_ids, str):
            rxn_ids = [rxn_ids]
         
        rxn2kcat = {}   
        if rxn_ids is None:
            # return all kcat values
            rxn2kcat = self.rxn2kcat
        else:
            # return reaction specific kcat values
            for rxn_id in rxn_ids:
                catalytic_event_id = self.catalytic_event_id.format(rxn_id)
                if catalytic_event_id not in self.catalytic_events:
                    warn('No reaction {0} found'.format(rxn_id))
                else:
                    # get kcat values
                    rxn2kcat[rxn_id] = self.rxn2kcat[rxn_id]
                    
        return rxn2kcat

    def remove_catalytic_event(self, catalytic_event: Union[CatalyticEvent, str]):
        """
        Function to remove a catalytic event from an enzyme
        Parameters
        ----------
        catalytic_event: CatalyticEvent or str
            catalytic event or identifier to remove
        """
        if isinstance(catalytic_event, str):
            try:
                catalytic_event = self.catalytic_events.get_by_id(catalytic_event)
            except:
                print(f'Catalytic event {catalytic_event} is not related to this enzyme and can thus not be removed!')

        #remove the event from the DictList
        self.catalytic_events.remove(catalytic_event)

        
class EnzymeComplex(Enzyme):
    """Upper level EnzymeComplex object containing information about the enzymes in a complex
       and link to the enzyme variables (CatalyticEvents) for each reaction the
       enzyme complex catalyzes.
       This class is used to generate enzyme instances from kcat values and contains
       information about the forward as well as the backward catalysis.

       There are two scenarios:
       - Promiscuous enzymes: a single enzyme complex can catalyze multiple reactions
       - Other: a single enzyme complex catalyzes a single reaction

       Parameters
       -------
       id : str
           Identifier for the enzyme complex (e.g. Uniprot ID)
       enzymes: DictList of cobra.core.Enzyme
            Enzyme objects associated with the enzyme complex
       rxn2kcat: Dict
           Dictionary with reaction ID, kcat value pairs for the forward (f)
           and backward (b) reaction
           (Example: {'PGI': {'f': 30, 'b': 0.1}})
       upper_bound: float
           Upper bound for the enzyme variable (default 1000.0)
       name: str
           Name of the enzyme (default None)
       molmass: float
           Molar mass of the enzyme (default 3.947778784340140e04)

       """

    # constant parameters
    DEFAULT_ENZYME_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]

    def __init__(
            self,
            id: str,
            enzymes: DictList,
            rxn2kcat: Dict,
            upper_bound: Union[int, float] = 1000.0,
            name: Optional[str] = None,
            molmass: Union[int, float] = DEFAULT_ENZYME_MOL_MASS,
    ):
        super().__init__(
            id = id,
            rxn2kcat=rxn2kcat,
            upper_bound = upper_bound,
            name=name,
            molmass=molmass
        )

        self.enzymes = None
        self.add_enzymes(enzymes)

    def add_enzymes(self, enzymes: DictList):
        for enzyme in enzymes:
            self.enzymes.append(enzyme)
            enzyme.enzyme_complex.append(self.id)
            self.molmass += enzyme.molmass