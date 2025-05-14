import math
from warnings import warn
from copy import copy, deepcopy

import pandas as pd
from cobra import Object
import os

from .configuration import Config
from .EnzymeSectors import EnzymeSector


class MembraneSector(EnzymeSector):
    def __init__(
            self,
            area_avail_0, #μm2
            area_avail_mu, #μm2/h
            alpha_numbers_dict: {},
            enzyme_location: {},
            cog_class: {} = None,
            max_area: float = 0.27,
            r_alpha: float = 0.00023, # radius of one alpha helix [um]
            cdw_per_cell: float = 0.28 * 1e-12,  # 0.28 pg
            n_a: float = 6.02214076 * 1e23,  # avogadro number
            configuration=Config):

        self.id = 'MembraneSector'
        self.alpha_numbers_dict = alpha_numbers_dict
        self.cog_class = cog_class
        self.enzyme_location = enzyme_location
        self.area_alpha = math.pi * math.pow((r_alpha), 2)  # area per alpha helix unit [um2]
        self.max_membrane_area = max_area #percentage of membrane area that can be covered by proteins
        self.unit_factor = 1e-3 * cdw_per_cell * n_a

        #Defining the slope and intercept
        self.intercept = area_avail_0 #μm2
        self.slope = area_avail_mu #μm2/h

    def add(self, model):

        print("Add membrane protein sector \n")
        model.membrane_sector = self
        self._add_membrane_constraint(model)
        pass

    def _add_membrane_constraint(self, model):

        self.membrane_proteins = {}

        coefficients = {
            model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable: -self.slope
        }

        for enz_complex in model.enzyme_variables:
            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)

            # Save membrane proteins (id, kcat and alpha number) inside of a dictionary
            if not isinstance(enz_complex, str) and alpha_number_for_complex != 0:
                self.membrane_proteins[enz_complex.id] = [enz_complex.kcats, alpha_number_for_complex]

            coeff = self._get_coeff_value(alpha_number_for_complex)

            coefficients[enz_complex.forward_variable] = coeff / self.max_membrane_area
            coefficients[enz_complex.reverse_variable] = coeff / self.max_membrane_area

        occupied_membrane = model.problem.Constraint(0, lb=0, ub=self.intercept, name='membrane')
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
            memprot_w_area_df = self.get_df_w_memprot_area(model, occupied_area)

            data_path = os.path.join('Results/PAM_parametrizer/Files/2025_03_11/memprot_data.xlsx')
            with pd.ExcelWriter(data_path, engine='openpyxl', mode='a') as writer:
                # Write the new DataFrame to a new sheet
                memprot_w_area_df.to_excel(writer, sheet_name=f'memprot_area_2', index=True)

        available_area = self.slope * model.objective.value + self.intercept

        return occupied_area, available_area

    def get_df_w_memprot_area(self, model, occupied_area):
        memprot_w_area = []


        for enz_complex in model.enzyme_variables:
            enz_complex_concentration = enz_complex.forward_variable.primal + enz_complex.reverse_variable.primal
            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)
            coeff = self._get_coeff_value(alpha_number_for_complex)

            for rxn_id, flux_dict in enz_complex.kcats.items():
                memprot_w_area.append({
                    'enzyme_id': enz_complex.id,
                    'Reaction': rxn_id,
                    'Forward Flux': flux_dict['f'] if 'f' in flux_dict else 0,
                    'Backward Flux': flux_dict['b'] if 'b' in flux_dict else 0,
                    'Occupied Area um2': coeff * enz_complex_concentration,
                    'Occupied Area %': coeff * enz_complex_concentration / occupied_area * 100,
                    'Contribution to protein pool': enz_complex_concentration * 1e-9 * enz_complex.molmass / model.p_tot * 100
                })

        memprot_w_area_df = pd.DataFrame(memprot_w_area)

        return memprot_w_area_df

    def change_available_membrane_area(self, new_max_area: float, model):
        self._update_membrane_constraint(new_max_area, model)

    def _update_membrane_constraint(self, new_max_area:float, model):
        self.max_membrane_area = new_max_area
        self.membrane_proteins = {}

        coefficients = {
            model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable: -self.slope
        }

        for enz_complex in model.enzyme_variables:
            alpha_number_for_complex = self._get_alpha_number_for_enz_complex(enz_complex)
            coeff = self._get_coeff_value(alpha_number_for_complex)

            coefficients[enz_complex.forward_variable] = coeff / new_max_area
            coefficients[enz_complex.reverse_variable] = coeff / new_max_area

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
            if enz in self.alpha_numbers_dict.keys() and self.enzyme_location[enz] == 'Cell membrane':
                alpha_numbers_in_complex.append(self.alpha_numbers_dict[enz])

        alpha_number_for_enz_complex = max(alpha_numbers_in_complex)

        return alpha_number_for_enz_complex

    def _get_coeff_value(self, alpha_number_for_complex:int):

        coeff = (1e-6 # correction for the solver issue
                 * alpha_number_for_complex
                 * self.area_alpha
                 * self.unit_factor)

        return coeff