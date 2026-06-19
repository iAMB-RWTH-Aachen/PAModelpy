import math
from warnings import warn
from copy import copy, deepcopy

import pandas as pd
from cobra import Object
import os
from typing import Union

from .configuration import Config
from .EnzymeSectors import EnzymeSector
from .Enzyme import Enzyme


class MembraneSector(EnzymeSector):
    """
    Membrane protein sector for constraining the occupancy of the
    inner membrane by integral membrane proteins (IMPs).

    The membrane occupancy constraint is formulated as:

        (S/V_available) * f >= Σ(
            E_conc
            * 10^3
            * CDW/V
            * N_A
            * N_alpha
            * A_alpha
        )

    where:

        - S/V_available : Available surface-to-volume ratio [µm²/fL]
        - f             : Fraction of membrane area available for proteins [-]
        - E_conc        : Enzyme concentration [mmol/gDW]
        - CDW/V         : Cell dry weight per volume [gDW/fL]
        - N_A           : Avogadro constant [1/mol]
        - N_alpha       : Number of transmembrane alpha helices [-]
        - A_alpha       : Area occupied by one alpha helix [µm²]

    The available membrane area is normalized by the volume and modeled as a linear function
    of the growth rate (µ):

        S/V_available = sv_0 + sv_slope * µ

    where:

        - sv_0      : Surface-to-volume ratio intercept [µm²/fL]
        - sv_slope  : Surface-to-volume ratio slope [µm²·h/fL]

    Args:
        sv_0 (float or int): Intercept of the surface-to-volume ratio relation [µm2/fL]
        sv_slope: (float or int): Slope of the surface-to-volume ratio relation [µm2.h/fL]
        alpha_numbers_dict (dict): Dictionary mapping enzyme IDs to the number of transmembrane alpha helices.
        enzyme_location (dict): Dictionary containing the cellular localization of enzymes.
        cog_class (dict): Dictionary containing COG classifications for enzymes.
        usable_area_fraction (float or int): Fraction of the total membrane area available for membrane proteins
        a_alpha (float or int): Membrane area occupied by a single alpha helix [µm2]
        cdw_per_volume (float or int): Cell dry weight per volume [gDW/fL]
        n_a (float or int): Avogadro constant [1/mol]
        separate_memprot_from_tpc (bool): If True, membrane proteins are excluded from the total protein constraint
        configuration (Config): PAModelpy configuration object

    """
    def __init__(
            self,
            sv_0: Union[int, float],
            sv_slope: Union[int, float], 
            alpha_numbers_dict: {},
            enzyme_location: {},
            cog_class: {} = None,
            usable_area_fraction: Union[int, float] = 0.5154,
            a_alpha: Union[int, float] = 1.4 * 1e-6, 
            cdw_per_volume: Union[int, float] = 268.36 * 1e-15,  
            n_a: Union[int, float] = 6.02214076 * 1e23,  
            separate_memprot_from_tpc: bool = True,
            unused_membrane_sector: object = None,
            configuration=Config):

        self.id = 'MembraneSector'
        self.alpha_numbers_dict = alpha_numbers_dict
        self.cog_class = cog_class
        self.enzyme_location = enzyme_location
        self.area_alpha = a_alpha
        self.usable_area_fraction = usable_area_fraction 
        self.unit_factor = 1e-3 * cdw_per_volume * n_a
        self.membrane_proteins = {}
        self._unused_membrane_sector = unused_membrane_sector

        # Defining the slope and intercept
        self.intercept = sv_0 #μm2
        self.slope = sv_slope #μm2/h

        self.separate_memprot_from_tpc = separate_memprot_from_tpc
    
    @property
    def unused_membrane_sector(self):
        return self._unused_membrane_sector
    
    @unused_membrane_sector.setter
    def unused_membrane_sector(self, unused_membrane_sector):
        self._unused_membrane_sector = unused_membrane_sector

        if unused_membrane_sector is not None:
            self.unused_membrane_sector._link_unused_enzyme_sector_to_membrane_sector()

    def add(self, model):

        print(f"Add membrane protein sector with {self.usable_area_fraction*100}% max inner membrane area\n")
        model.membrane_sector = self
        self._add_membrane_constraint(model)

        # Add the metabolic model to unused_membrane_sector object if not None
        if self.unused_membrane_sector is not None:
            self.unused_membrane_sector.model = model
            self.unused_membrane_sector._link_unused_enzyme_sector_to_membrane_sector(model=model)
        pass

    def _add_membrane_constraint(self, model):

        coefficients = {
            model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable: -self.slope * self.usable_area_fraction
        }

        for enz_complex in model.enzyme_variables:

            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)

            # Save membrane proteins (id, kcat and alpha number) inside of a dictionary
            if not isinstance(enz_complex, str) and alpha_number_for_complex != 0:
                self.membrane_proteins[enz_complex.id] = [enz_complex.kcats, alpha_number_for_complex]

            coeff = self._get_coeff_value(alpha_number_for_complex)

            coefficients[enz_complex.forward_variable] = coeff
            coefficients[enz_complex.reverse_variable] = coeff 

            # Exclude membrane enzymes from the total protein constraint
            if self.separate_memprot_from_tpc and enz_complex.id in self.membrane_proteins:
                model.exclude_variable_from_constraint(variable=enz_complex, constraint_name=model.TOTAL_PROTEIN_CONSTRAINT_ID)

        occupied_membrane = model.problem.Constraint(0, lb=0, ub=self.intercept*self.usable_area_fraction, name='membrane')
        model.add_cons_vars(occupied_membrane)
        model.solver.update()
        occupied_membrane.set_linear_coefficients(coefficients=coefficients)

    def calculate_occupied_membrane(self,
                                    model,
                                    get_df:bool = False,
                                    ):
        occupied_area = 0
        memprot_conc = 0

        # Occupied area from membrane proteins
        for enz_complex in model.enzyme_variables:
            enz_complex_concentration = (
                enz_complex.forward_variable.primal
                + enz_complex.reverse_variable.primal
            )

            alpha_number_for_complex = (
                self._get_alpha_number_for_enz_complex(enz_complex)
            )

            coeff = self._get_coeff_value(alpha_number_for_complex)

            occupied_area += coeff * enz_complex_concentration

        # Add occupancy from the unused membrane sector if enabled
        unused_area = 0

        if self.unused_membrane_sector is not None:
            ups_sector = model.sectors.get_by_id('UnusedEnzymeSector')
            lin_rxn = model.reactions.get_by_id(ups_sector.id_list[0])

            unused_flux = (
                lin_rxn.forward_variable.primal
                - lin_rxn.reverse_variable.primal
            )

            unused_area = (
                self.unused_membrane_sector.intercept
                + self._unused_membrane_sector.slope * unused_flux
            )

        total_occupied_area = occupied_area + unused_area

        if get_df:
            df = self.get_prot_occupancy_df(model, total_occupied_area)
            return df

        # Total membrane area available from geometry
        total_available_area = (
            self.slope * model.objective.value
            + self.intercept
        ) * self.usable_area_fraction

        return total_occupied_area, total_available_area

   

    def get_prot_occupancy_df(self, model, occupied_area):
        df = []

        for enz_complex in model.enzyme_variables:
            enz_complex_concentration = enz_complex.forward_variable.primal + enz_complex.reverse_variable.primal
            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)
            coeff = self._get_coeff_value(alpha_number_for_complex)

            for rxn_id, flux_dict in enz_complex.kcats.items():
                df.append({
                    'enzyme_id': enz_complex.id,
                    'Reaction': rxn_id,
                    'Forward Flux': flux_dict['f'] if 'f' in flux_dict else 0,
                    'Backward Flux': flux_dict['b'] if 'b' in flux_dict else 0,
                    'Alpha Number': alpha_number_for_complex,
                    'Occupied Area um2': coeff * enz_complex_concentration,
                    'Occupied Area %': coeff * enz_complex_concentration / occupied_area * 100,
                    'Contribution to protein pool': enz_complex_concentration * 1e-9 * enz_complex.molmass / model.p_tot * 100
                })

        df = pd.DataFrame(df)

        return df

    def change_available_membrane_area(self, new_max_area: float, model):
        self._update_membrane_constraint(new_max_area, model)

    def _update_membrane_constraint(self,
                                    new_usable_area_fraction: float,
                                    model
                                    ):
        """
        Fast update of membrane constraint:
        Only updates terms that depend on usable membrane area fraction.

        Assumes enzyme coefficients are already initialized and do NOT change here.
        """

        self.usable_area_fraction = new_usable_area_fraction
        membrane_constraint = model.constraints["membrane"]

        # biomass coefficient update
        biomass_var = model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable
        membrane_constraint.set_linear_coefficients({biomass_var: -self.slope * new_usable_area_fraction})

        # upper bound (intercept term)
        ub = self.intercept * new_usable_area_fraction
        membrane_constraint.ub = ub

        # unused membrane sector (if enabled)
        if self.unused_membrane_sector is not None:
            self.unused_membrane_sector._link_unused_enzyme_sector_to_membrane_sector(model)

        model.solver.update()

        return self

    def _get_alpha_number_for_enz_complex(self, enz_complex):
        if isinstance(enz_complex, str):
            enzymes = enz_complex.split("_")
        else:
            enzymes = enz_complex.id.split("_")
        alpha_numbers_in_complex = [0]  # zero if enzyme is not in membrane

        for enz in enzymes:
            if enz in self.alpha_numbers_dict.keys() and self.enzyme_location[enz] == 'Cell inner membrane':
                alpha_numbers_in_complex.append(self.alpha_numbers_dict[enz])

        alpha_number_for_enz_complex = sum(alpha_numbers_in_complex)

        return alpha_number_for_enz_complex

    def _get_coeff_value(self, alpha_number_for_complex:int):
        '''
        Returns coefficient value in um2*gDW/mmol*fL
        '''

        coeff = (1e-6 # correction for the solver issue
                 * alpha_number_for_complex
                 * self.area_alpha
                 * self.unit_factor)

        return coeff

class UnusedMembraneSector(EnzymeSector):
    DEFAULT_ALPHA_NUMBER = 12  # default alpha helix unit number for transport proteins [-] (PJ Henderson 1993)
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]
    DEFAULT_TOTAL_PROTEIN_CONCENTRATION = 0.258 # global protein concentration of an E.coli cell 

    # class with all the information on 'excess' membrane sector (linear dependent on substrate uptake rate)
    def __init__(self, model=None):
        self.mol_mass = [self.DEFAULT_ALPHA_NUMBER]
        self._model = model
        self.intercept = 0
        self.slope = 0

        if model is not None:
            ups_slope, ups_intercept = self._get_ue_sector_parameters(model)
            self.set_sector_slope_and_intercept_from_ups_sector(model, ups_slope, ups_intercept)


    @property
    def model(self):
        return self._model
    
    @model.setter
    def model(self, model):
        self._model = model

        ups_slope, ups_intercept = self._get_ue_sector_parameters(model) 

        # Convert the slope and intercept to the correct unit 
        ups_slope, ups_intercept = self._get_ue_sector_parameters(model)
        self.set_sector_slope_and_intercept_from_ups_sector(model, ups_slope, ups_intercept)

    def set_sector_slope_and_intercept_from_ups_sector(self, model, ups_slope, ups_intercept):
        """
        Take the unused protein sectors's (ups) intercept and slope, convert it to suitable unit for unused membrane sector, and finally update unused membrane sector's 
        intercept and slope
        """
        conversion_unit = self._get_conversion_unit()
        unused_membrane_fraction = 0.15
        # unused_membrane_fraction = (self.DEFAULT_TOTAL_PROTEIN_CONCENTRATION - model.p_tot) / self.DEFAULT_TOTAL_PROTEIN_CONCENTRATION # fraction of unused enzyme that can be allocated to the membrane 
        self.intercept = ups_intercept * conversion_unit * unused_membrane_fraction
        self.slope = ups_slope * conversion_unit *unused_membrane_fraction

    def _get_ue_sector_parameters(self, model):
        """
        Get UnusedEnzymeSector's slope and intercept for the given model

        Args:
            model(Model, PAModel): Constraint-based metabolic model

        Returns:
            ups_intercept = Unused protein sector's intercept (Amount of protein allocated to the excess enzyme sector at zero substrate uptake in g/gDW)
            ups_slope = Unused protein sector's slope (Slope of linear relation with growth/substrate uptake in g/gDW/h)
        """
        if isinstance(model.sectors.get_by_id('UnusedEnzymeSector').ups_0, list):
            ups_intercept =  model.sectors.get_by_id('UnusedEnzymeSector').ups_0[0]
        else:
            self.ups_intercept = model.sectors.get_by_id('UnusedEnzymeSector').ups_0# amount of protein allocated to the excess enzyme sector at zero substrate uptake in g/gDW
        ups_slope = model.sectors.get_by_id('UnusedEnzymeSector').ups_mu # slope of linear relation with growth/substrate uptake in g/gDW/h

        return ups_slope, ups_intercept

    def _get_conversion_unit(self):
        """
        Returns conversion unit in um2*gDW/(g*fL) to transform unused protein sector's slope/intercept to unused membrane sector slope/intercept
        """
        # Get initial coefficient value from the existing method in MembraneSector with unit um2*gDW/mmol*fL
        coeff_value = self.model.sectors.get_by_id('MembraneSector')._get_coeff_value(self.DEFAULT_ALPHA_NUMBER)

        # Transfrom the coefficient value to conversion unit of um2*gDW/g*fL
        conversion_unit = coeff_value / self.DEFAULT_MOL_MASS * 1e9 # multiply with 1e9 to correct enzyme unit and solver tolerance

        return conversion_unit
    
    def _link_unused_enzyme_sector_to_membrane_sector(self, model=None):
        model = model if model is not None else self.model
        ups_sector = model.sectors.get_by_id('UnusedEnzymeSector') # unused protein sector
        lin_rxn = model.reactions.get_by_id(ups_sector.id_list[0])
        
        # relate the sector to the total membrane sector if it's there
        membrane_constraint = model.constraints['membrane']
        if 'membrane' in model.constraints.keys():
            # add parts of constraint corresponding to the unused enzyme sector to the total_membrane constraint
            # 1. subtract the intercept value from the sum of surface to volume ratio in the membrane (S/V_available == S/V_total - S/V_unused)
            membrane_ub = membrane_constraint.ub
            membrane_constraint.ub = (
                membrane_ub - self.intercept
            )

            # 2. add the slope and variable of unused membrane sector to the right hand side of the membrane constraint
            for direction_variable, coeff in zip([lin_rxn.forward_variable, lin_rxn.reverse_variable], [self.slope, -self.slope]):

                membrane_constraint.set_linear_coefficients(
                    {
                        direction_variable: coeff
                    })