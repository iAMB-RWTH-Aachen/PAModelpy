import math
from warnings import warn
from copy import copy, deepcopy

import pandas as pd
from cobra import Object
import os

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
            configuration=Config):

        self.id = 'MembraneSector'
        self.alpha_numbers_dict = alpha_numbers_dict
        self.cog_class = cog_class
        self.enzyme_location = enzyme_location
        self.area_alpha = a_alpha
        self.usable_area_fraction = usable_area_fraction 
        self.unit_factor = 1e-3 * cdw_per_volume * n_a
        self.membrane_proteins = {}

        #Defining the slope and intercept
        self.intercept = sv_0 #μm2
        self.slope = sv_slope #μm2/h

        self.separate_memprot_from_tpc = separate_memprot_from_tpc

    def add(self, model):

        print(f"Add membrane protein sector with {self.usable_area_fraction*100}% max inner membrane area\n")
        model.membrane_sector = self
        self._add_membrane_constraint(model)
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

    def calculate_occupied_membrane(self, model, get_df:bool = False, get_memprot_contribution:bool = False):
        occupied_area = 0
        memprot_conc = 0

        for enz_complex in model.enzyme_variables:
            enz_complex_concentration = enz_complex.forward_variable.primal + enz_complex.reverse_variable.primal
            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)
            coeff = self._get_coeff_value(alpha_number_for_complex)
            occupied_area += coeff * enz_complex_concentration

            if alpha_number_for_complex != 0:
                memprot_conc += enz_complex_concentration * 1e-9 * enz_complex.molmass # g_enz/g_DW

        if get_memprot_contribution: # Membrane protein contribution to the total protein pool
            memprot_contribution = memprot_conc / model.p_tot * 100 # Result in %

            return memprot_contribution

        if get_df:
            df = self.get_prot_occupancy_df(model, occupied_area)

            return df

        available_area = (self.slope * model.objective.value + self.intercept) * self.usable_area_fraction # available membrane area for inner membrane proteins

        return occupied_area, available_area

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

    def _update_membrane_constraint(self, new_max_area:float, model):
        self.usable_area_fraction = new_max_area
        self.membrane_proteins = {}

        coefficients = {
            model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable: -self.slope * new_max_area
        }

        for enz_complex in model.enzyme_variables:
            if self.separate_memprot_from_tpc and enz_complex.id in self.membrane_proteins:
                var_names = {v.name for v in model.constraints[model.TOTAL_PROTEIN_CONSTRAINT_ID].expression.free_symbols}
                if any(name.startswith(enz_complex.id) for name in var_names):
                    model.exclude_variable_from_constraint(variable=enz_complex, constraint_name=model.TOTAL_PROTEIN_CONSTRAINT_ID)

            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)
            coeff = self._get_coeff_value(alpha_number_for_complex)

            coefficients[enz_complex.forward_variable] = coeff 
            coefficients[enz_complex.reverse_variable] = coeff 

        model.constraints['membrane'].ub = self.intercept*new_max_area
        model.constraints['membrane'].set_linear_coefficients(coefficients=coefficients)
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

        coeff = (1e-6 # correction for the solver issue
                 * alpha_number_for_complex
                 * self.area_alpha
                 * self.unit_factor)

        return coeff
