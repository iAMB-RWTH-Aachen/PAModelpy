import pandas as pd
from tabulate import tabulate
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import cobra
import sys
import numpy as np

from IPython.display import display
from IPython import get_ipython

from .configuration import Config

from typing import Union


class PAMValidator(object):
    RESULT_DIR = os.path.join(os.path.split(os.getcwd())[0], "Results")
    OUT_FILE = os.path.join(RESULT_DIR, "flux_rates_vs_glc.xls")
    GLUCOSE_EXCHANGE_RXNID = Config.GLUCOSE_EXCHANGE_RXNID
    OXYGEN_UPTAKE_RXNID = Config.OXYGEN_UPTAKE_RXNID
    BIOMASS_REACTION = Config.BIOMASS_REACTION
    ACETATE_EXCRETION_RXNID = Config.ACETATE_EXCRETION_RXNID
    CO2_EXHANGE_RXNID = Config.CO2_EXHANGE_RXNID
    PHYS_RXN_IDS = Config.PHYS_RXN_IDS
    MW_GLC = 180.15588  # g/mol
    MW_ACETATE = 59.04
    MW_CO2 = 44.009
    GRADIENT_MAX = 11  # mmol/gdw/h
    GRADIENT_STEP = 0.5  # mmol/gdw/h
    GRADIENT_MIN = 0  # mmol/gdw/h

    def __init__(self, model, phys_file, configuration=Config):
        self.GLUCOSE_EXCHANGE_RXNID = configuration.GLUCOSE_EXCHANGE_RXNID
        self.OXYGEN_UPTAKE_RXNID = configuration.OXYGEN_UPTAKE_RXNID
        self.BIOMASS_REACTION = configuration.BIOMASS_REACTION
        self.ACETATE_EXCRETION_RXNID = configuration.ACETATE_EXCRETION_RXNID
        self.CO2_EXHANGE_RXNID = configuration.CO2_EXHANGE_RXNID
        self.PHYS_RXN_IDS = configuration.PHYS_RXN_IDS

        self.model = model
        self.physiology_data = pd.read_excel(phys_file, sheet_name="Fluxes")
        self.yield_data = pd.read_excel(phys_file, sheet_name="Yields")

        self.results = pd.DataFrame(
            columns=["Y_biomass", "Y_acetate", "Y_co2"] + self.PHYS_RXN_IDS
        )
        self._additional_rxns = []

    @staticmethod
    def check_kernel_type():
        try:
            ipy_str = str(type(get_ipython()))
            if "zmqshell" in ipy_str:
                # jupyter
                return True
            if "terminal" in ipy_str:
                # ipython
                return False
        except:
            # terminal
            return False

    def add_rxn_to_save(self, rxn_list: Union[str, list]):
        if isinstance(rxn_list, str):
            rxn_list = [rxn_list]

        self.PHYS_RXN_IDS += rxn_list
        for rxn in rxn_list:
            self.results[rxn] = 0
        self._additional_rxns.append(rxn_list)

    def validate_literature(self):
        self.parse_data()
        self.run_simulations()
        # for pretty printing: check if the code is ran from a jupyter notebook or not
        jupyter = self.check_kernel_type()

        if jupyter:
            pd.set_option("display.max_rows", None)
            display(self.results)
        else:
            print(tabulate(self.results, headers="keys", tablefmt="psql"))
        self.plot_mu()

    def validate_range(
        self,
        c_uptake_rxn: str = GLUCOSE_EXCHANGE_RXNID,
        pfba: bool = False,
        show: bool = True,
        save: bool = False,
    ):
        # reset result dataframe
        self.results = pd.DataFrame(
            columns=["Y_biomass", "Y_acetate", "Y_co2"] + self.PHYS_RXN_IDS
        )
        try:
            self.parse_data()
        except:
            pass
        self.run_simulations_gradient(c_uptake_rxn=c_uptake_rxn, pfba=pfba)
        jupyter = self.check_kernel_type()

        if jupyter:
            pd.set_option("display.max_rows", None)
            display(self.results)
        else:
            print(tabulate(self.results, headers="keys", tablefmt="psql"))
        if show:
            self.plot_vs_mu(c_uptake_rxn)
            if self._additional_rxns != []:
                self.custom_plot(self._additional_rxns, c_uptake_rxn=c_uptake_rxn)
        if save:
            # make result directory if it not already exists
            if not os.path.exists(self.RESULT_DIR):
                os.mkdir(self.RESULT_DIR)
            self.results.to_excel(self.OUT_FILE, engine="openpyxl")

    def parse_data(self):
        self.physiology_data = pd.melt(
            self.physiology_data,
            id_vars="Reaction_ID",
            var_name="Reference",
            value_vars=[
                "Rijsewijk_2011",
                "Nanchen_2006_1",
                "Nanchen_2006_2",
                "Nanchen_2006_3",
                "Nanchen_2006_4",
            ],
        )

    def run_simulations(self):
        i = 0
        print("Running simulations to validate the model...\n")
        # physiology dataset
        for ref in pd.unique(self.physiology_data.Reference):
            # set glucose uptake rate
            df = self.physiology_data[self.physiology_data.Reference == ref]
            bound = df[df.Reaction_ID == self.GLUCOSE_EXCHANGE_RXNID].value.iloc[0]
            self.set_glc_uptake_bounds(bound)
            self.model.optimize()
            if self.model.solver.status == "optimal":
                self.store_results()
        print("Still validating...\n")
        for ref in pd.unique(self.yield_data.Reference):
            df = self.yield_data[self.yield_data.Reference == ref]
            for glc_flux in pd.unique(df.Glucose_IN):
                self.set_glc_uptake_bounds(glc_flux)
                self.model.optimize()
                if self.model.solver.status == "optimal":
                    self.store_results()
                if i == 2:
                    cobra.core.solution.get_solution(self.model).fluxes.to_csv(
                        os.path.join(
                            os.path.split(os.getcwd())[0], "Results", "test_out_mu.tsv"
                        ),
                        sep="\t",
                    )
                i += 1

    def run_simulations_gradient(
        self, pfba: bool = False, c_uptake_rxn: str = GLUCOSE_EXCHANGE_RXNID
    ):
        # run simulations for increasing glucose consumption rates
        print("Running simulations to validate the model...\n")
        self.model.reset_objective()
        for flux in np.arange(self.GRADIENT_MIN, self.GRADIENT_MAX, self.GRADIENT_STEP):
            self.model.change_reaction_bounds(c_uptake_rxn, -flux, -flux)
            # if glucose is not the substrate turn off glucose uptake reaction
            if c_uptake_rxn != self.GLUCOSE_EXCHANGE_RXNID:
                self.set_glc_uptake_bounds(0)
            if pfba:
                self.model.reset_objective()
            self.model.optimize()
            if self.model.solver.status == "optimal":
                if pfba:
                    self.model.pfba(
                        fraction_of_optimum=0.1, proteins=True, reactions=False
                    )
                self.store_results(c_uptake_rxn=c_uptake_rxn)
            print("Done running simulations\n")

    def run_simulations_glc_o2_gradient(
        self,
        oxygen_gradient: list,
        params_to_save: Union[str, list] = "R_TranslationalProteinSector",
    ):
        """
        Function to run simulations of different oxygen gradients for a range of growth rates.

        This will simulate growth for the entire range of glucose concentrations for each oxygen uptake rate as given by the input.

        Parameters:
            oxygen_gradient (list): List of upper bounds for the oxygen uptake reaction to loop over.
            params_to_save (optional): string or list, which parameter(s) to save for further analysis (default: translational protein sector constraint).

        Returns:
            results (list of dataframes): Saves the growth rate, glucose uptake rate, and the user-defined parameters for each oxygen uptake rate in separate dataframes.
        """

        print("Running simulations of glucose and oxygen gradients...\n")
        result = list()
        if isinstance(params_to_save, str):
            params_to_save = [params_to_save]
        rxns_to_save = [
            self.BIOMASS_REACTION,
            self.GLUCOSE_EXCHANGE_RXNID,
            self.OXYGEN_UPTAKE_RXNID,
        ] + params_to_save
        for o2 in oxygen_gradient:
            df = pd.DataFrame(columns=rxns_to_save)
            self.set_rxn_uptake_bounds(self.OXYGEN_UPTAKE_RXNID, o2)
            print(self.model.reactions.get_by_id(self.OXYGEN_UPTAKE_RXNID).bounds)
            for glc_flux in range(0, self.GRADIENT_MAX, self.GRADIENT_STEP):
                self.set_glc_uptake_bounds(glc_flux / 10)
                self.model.reset_objective()
                self.model.optimize()
                if self.model.solver.status == "optimal":
                    df = self.store_manual_results(rxns_to_save, df)
            result += [df]
        return result

    def run_simulations_ups(
        self,
        ups_gradient: list,
        params_to_save: Union[str, list] = "R_TranslationalProteinSector",
    ):
        """
        Function to run simulations with increasing unused enzyme sectors proportions for a range of growth rates.

        This will simulate growth for the entire range of glucose concentrations for a range of fractions of ups_0 as given by the input.

        Parameters:
            ups_gradient (list): List of upper bounds for the oxygen uptake reaction to loop over.
            params_to_save (optional): string or list, which parameter(s) to save for further analysis (default: translational protein sector constraint).

        Returns:
            results (list of dataframes): Saves the growth rate, glucose uptake rate, and the user-defined parameters for each oxygen uptake rate in separate dataframes.
        """

        print(
            "Running simulations of different unused protein intercepts and glucose gradients...\n"
        )
        result = list()
        if isinstance(params_to_save, str):
            params_to_save = [params_to_save]
        rxns_to_save = [
            self.BIOMASS_REACTION,
            self.GLUCOSE_EXCHANGE_RXNID,
            self.OXYGEN_UPTAKE_RXNID,
        ] + params_to_save
        unused_protein_sector = self.model.sectors.get_by_id("UnusedProteinSector")
        ups_0 = unused_protein_sector.ups_0[0]
        for perc in ups_gradient:
            df = pd.DataFrame(columns=rxns_to_save)
            unused_protein_sector.ups_0[0] = ups_0 * perc
            # self.set_rxn_uptake_bounds(self.OXYGEN_UPTAKE_RXNID, perc)
            for glc_flux in range(0, self.GRADIENT_MAX, self.GRADIENT_STEP):
                self.set_glc_uptake_bounds(glc_flux / 10)
                self.model.reset_objective()
                self.model.optimize()
                if self.model.solver.status == "optimal":
                    df = self.store_manual_results(rxns_to_save, df)
            result += [df]
        return result

    def set_glc_uptake_bounds(self, bound: int):
        rxn_glc_in = self.model.reactions.get_by_id(self.GLUCOSE_EXCHANGE_RXNID)
        # check if model is reversible
        if self.GLUCOSE_EXCHANGE_RXNID[-1] == "b":
            rxn_glc_in._lower_bound, rxn_glc_in._upper_bound = bound, bound
        else:
            rxn_glc_in._lower_bound, rxn_glc_in._upper_bound = -bound, -bound
        rxn_glc_in.update_variable_bounds()

    def set_rxn_uptake_bounds(
        self, rxn_id: str, bound: int, upper: bool = True, lower: bool = False
    ):
        rxn = self.model.reactions.get_by_id(rxn_id)
        if rxn_id[-1] == "b":
            if upper:
                rxn._upper_bound = bound
            if lower:
                rxn._lower_bound = bound
        else:
            if upper:
                rxn._lower_bound = -bound
            if lower:
                rxn._upper_bound = -bound
        rxn.update_variable_bounds()

    def store_results(self, c_uptake_rxn: str = GLUCOSE_EXCHANGE_RXNID):
        result = {}
        c_upt_flux = self.model.reactions.get_by_id(c_uptake_rxn).flux
        # store the relevant flux values
        if c_uptake_rxn != self.GLUCOSE_EXCHANGE_RXNID:
            for i, rxn in enumerate(self.PHYS_RXN_IDS):
                if rxn == self.GLUCOSE_EXCHANGE_RXNID:
                    self.PHYS_RXN_IDS[i] = c_uptake_rxn
                    self.results = pd.DataFrame(columns=self.PHYS_RXN_IDS)
                    # self.results = pd.DataFrame(columns=[ 'Y_biomass', 'Y_acetate', 'Y_co2']+self.PHYS_RXN_IDS)
        for id in self.PHYS_RXN_IDS:
            # getting reaction or variable
            try:
                result[id] = self.model.reactions.get_by_id(id).flux
            except:
                result[id] = self.model.variables[id].primal
            # calculate yields based on glucose for exchange reactions
            if (
                ("EX" in id or "BIOMASS" in id)
                and "glc" not in id
                and "_o2" not in id
                and abs(c_upt_flux) > 0
            ):
                if "BIOMASS" in id:
                    id_ext = "biomass"
                    MW = 1
                if "ac" in id:
                    id_ext = "acetate"
                    MW = 1  # self.MW_ACETATE * 1e-3 #[mmol/mol]
                if "_co2" in id:
                    id_ext = "co2"
                    MW = 1  # self.MW_CO2 * 1e-3 #[mmol/mol]
                yield_id = "Y_" + id_ext
                result[yield_id] = (self.model.reactions.get_by_id(id).flux * MW) / (
                    c_upt_flux * self.MW_GLC * 1e-3
                )
        # appending the results to the dataframe
        result_series = pd.Series(result)
        # sort the series to match the format of the output dataframe
        result_series_sorted = result_series[self.results.columns.to_list()]
        # add the resulting row to the result dataframe
        self.results.loc[len(self.results.index)] = result_series_sorted.to_list()

    def store_manual_results(self, vars_to_save: list, df: pd.DataFrame):
        row = dict()
        for id in vars_to_save:
            try:
                # for variables
                row[id] = self.model.variables[id].primal
            except:
                # for reactions:
                row[id] = self.model.reactions.get_by_id(id).flux
        # appending the results to the dataframe
        df = pd.concat(
            [df, pd.DataFrame({key: pd.Series(value) for key, value in row.items()})],
            ignore_index=True,
        )
        return df

    def plot_mu(self):
        print("Plotting the results against reference values. \n")
        # Initialize figure with subplots
        fig = make_subplots(
            rows=2,
            cols=2,
            specs=[
                [{"type": "scatter"}, {"type": "scatter"}],
                [{"type": "scatter"}, {"type": "scatter"}],
            ],
            subplot_titles=(
                "Growth rate",
                "Biomass yield",
                "CO2 excretion",
                "CO2 yield",
                "Acetate excretion",
                "Acetate yield",
            ),
        )
        result = self.results.__deepcopy__()
        result = result.merge(
            self.yield_data, left_on=self.GLUCOSE_EXCHANGE_RXNID, right_on="Glucose_IN"
        )
        # print(tabulate(result, headers='keys', tablefmt='psql'))
        # growth rate
        fig.add_trace(
            go.Scatter(
                x=result[self.BIOMASS_REACTION],
                y=result["mu"],
                mode="markers",
                showlegend=False,
                marker=dict(color="blue"),
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=[0, 0.5, 1.4],
                y=[0, 0.5, 1.4],
                mode="lines",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=1,
            col=1,
        )
        fig.update_xaxes(title_text="Simulated growth rate [1/h]", row=1, col=1)
        fig.update_yaxes(
            title_text="Experimentally determined growth rate [1/h]", row=1, col=1
        )

        # biomass yield
        fig.add_trace(
            go.Scatter(
                x=result["Y_biomass"],
                y=result["Yield_biomass"],
                mode="markers",
                showlegend=False,
                marker=dict(color="blue"),
            ),
            row=1,
            col=2,
        )
        fig.add_trace(
            go.Scatter(
                x=[0, 0.5, 0.6],
                y=[0, 0.5, 0.6],
                mode="lines",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=1,
            col=2,
        )
        fig.update_xaxes(title_text="Simulated biomass yield [g_x/g_glc]", row=1, col=2)
        fig.update_yaxes(
            title_text="Experimentally determined biomass yield [g_x/g_glc]",
            row=1,
            col=2,
        )

        # CO2 yield
        fig.add_trace(
            go.Scatter(
                x=result["Y_co2"],
                y=result["Yield_co2"],
                mode="markers",
                showlegend=False,
                marker=dict(color="blue"),
            ),
            row=2,
            col=2,
        )
        fig.add_trace(
            go.Scatter(
                x=[0, 15, 30],
                y=[0, 15, 30],
                mode="lines",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=2,
            col=2,
        )
        fig.update_xaxes(title_text="Simulated CO2 yield [mmol_x/g_glc]", row=2, col=2)
        fig.update_yaxes(
            title_text="Experimentally determined CO2 yield [mmol_co2/g_glc]",
            row=2,
            col=2,
        )

        # CO2 out
        fig.add_trace(
            go.Scatter(
                x=result[self.CO2_EXHANGE_RXNID],
                y=result["CO2_OUT"],
                mode="markers",
                showlegend=False,
                marker=dict(color="blue"),
            ),
            row=2,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=[0, 15, 30],
                y=[0, 15, 30],
                mode="lines",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=2,
            col=1,
        )
        fig.update_xaxes(
            title_text="Simulated CO2 exchange rate [mmol/g_cdw/h]", row=2, col=1
        )
        fig.update_yaxes(
            title_text="Experimentally determined CO2 exchange rate [mmol/g_cdw/h]",
            row=2,
            col=1,
        )
        # fig.write_image('Results/Images/mu_co2_plots.png') #depends on the kaleido or orca package

        fig.show()
        # fig.write_image('here.png')

    def plot_vs_mu(self, c_uptake_rxn: str = GLUCOSE_EXCHANGE_RXNID):
        print("Plotting the results against reference values. \n")
        # Initialize figure with subplots
        fig = make_subplots(
            rows=2,
            cols=2,
            specs=[
                [{"type": "scatter"}, {"type": "scatter"}],
                [{"type": "scatter"}, {"type": "scatter"}],
            ],
            subplot_titles=(
                "Glucose consumption",
                "Acetate production",
                "O2 consumption",
                "CO2 production",
            ),
        )

        self.results[c_uptake_rxn] = abs(self.results[c_uptake_rxn])
        self.results[self.OXYGEN_UPTAKE_RXNID] = abs(
            self.results[self.OXYGEN_UPTAKE_RXNID]
        )
        self.yield_data[c_uptake_rxn] = abs(self.yield_data[c_uptake_rxn])
        self.yield_data[self.OXYGEN_UPTAKE_RXNID] = abs(
            self.yield_data[self.OXYGEN_UPTAKE_RXNID]
        )
        self.yield_data[self.BIOMASS_REACTION] = abs(
            self.yield_data[self.BIOMASS_REACTION]
        )

        # glucose consumption
        fig.add_trace(
            go.Scatter(
                x=self.results[c_uptake_rxn],
                y=self.results[self.BIOMASS_REACTION],
                showlegend=False,
                mode="lines",
                marker=dict(color="blue"),
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=self.yield_data[c_uptake_rxn],
                y=self.yield_data[self.BIOMASS_REACTION],
                mode="markers",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=1,
            col=1,
        )
        fig.update_yaxes(title_text="Growth rate [1/h]", row=1, col=1)
        fig.update_xaxes(title_text="Glucose consumption [mmol/gdw/h]", row=1, col=1)

        # acetate excretion
        fig.add_trace(
            go.Scatter(
                x=self.results[c_uptake_rxn],
                y=self.results[self.ACETATE_EXCRETION_RXNID],
                showlegend=False,
                mode="lines",
                marker=dict(color="blue"),
            ),
            row=1,
            col=2,
        )
        fig.add_trace(
            go.Scatter(
                x=self.yield_data[c_uptake_rxn],
                y=self.yield_data[self.ACETATE_EXCRETION_RXNID],
                mode="markers",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=1,
            col=2,
        )
        fig.update_xaxes(title_text="Glucose consumption [mmol/gdw/h]", row=1, col=2)
        fig.update_yaxes(
            title_text="Acetate production rate [mmol/gdw/h]", row=1, col=2
        )

        # O2 uptake
        fig.add_trace(
            go.Scatter(
                x=self.results[c_uptake_rxn],
                y=self.results[self.OXYGEN_UPTAKE_RXNID],
                mode="lines",
                showlegend=False,
                marker=dict(color="blue"),
            ),
            row=2,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=self.yield_data[c_uptake_rxn],
                y=self.yield_data[self.OXYGEN_UPTAKE_RXNID],
                mode="markers",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=2,
            col=1,
        )
        fig.update_xaxes(title_text="Glucose consumption [mmol/gdw/h]", row=2, col=1)
        fig.update_yaxes(title_text="O2 consumption rate [mmol/gdw/h]", row=2, col=1)

        # CO2 excretion
        fig.add_trace(
            go.Scatter(
                x=self.results[c_uptake_rxn],
                y=self.results[self.CO2_EXHANGE_RXNID],
                mode="lines",
                showlegend=False,
                marker=dict(color="blue"),
            ),
            row=2,
            col=2,
        )
        fig.add_trace(
            go.Scatter(
                x=self.yield_data[c_uptake_rxn],
                y=self.yield_data[self.CO2_EXHANGE_RXNID],
                mode="markers",
                showlegend=False,
                line=dict(color="black", dash="dash", width=1),
            ),
            row=2,
            col=2,
        )
        fig.update_xaxes(title_text="Glucose consumption [mmol/gdw/h]", row=2, col=2)
        fig.update_yaxes(title_text="CO2 production rate [mmol/gdw/h]", row=2, col=2)

        fig.show()

    def custom_plot(
        self,
        rxn_ids: list,
        valid_dataframe: pd.DataFrame = None,
        xaxis: str = None,
        c_uptake_rxn: str = GLUCOSE_EXCHANGE_RXNID,
    ):
        """
        Function to plot the results of custom reactions.

        Parameters:
            rxn_ids (list of str): Reaction identifiers of the reactions to be plotted.
            valid_dataframe (pandas.DataFrame, optional): A DataFrame with experimental data to validate the results with.
                The columns should be the same as the rxn_id of the reaction to be plotted and the reaction which should be plotted
                on the x-axis (by default the glucose exchange reaction `EX_glc__D_e_b`). If the DataFrame is not provided,
                only the simulation results will be plotted.
            xaxis (str, optional): The reaction identifier of the reaction which should be plotted on the x-axis (default: `EX_glc__D_e_b`).

        Returns:
            Prints scatter plots of the model simulations vs. experimental data points (if provided).
        """

        if xaxis is None:
            xaxis = c_uptake_rxn
        fig = make_subplots(
            rows=1,
            cols=len(rxn_ids),
            specs=[[{"type": "scatter"}] * len(rxn_ids)],
            subplot_titles=tuple(rxn_id[0] for rxn_id in rxn_ids),
        )
        for i, rxn in enumerate(rxn_ids):
            fig.add_trace(
                go.Scatter(
                    x=self.results[xaxis],
                    y=self.results[rxn[0]],
                    mode="lines",
                    showlegend=False,
                    marker=dict(color="blue"),
                ),
                row=1,
                col=i + 1,
            )
            if valid_dataframe is not None:
                fig.add_trace(
                    go.Scatter(
                        x=valid_dataframe[xaxis],
                        y=valid_dataframe[rxn],
                        mode="markers",
                        showlegend=False,
                        line=dict(color="black", dash="dash", width=1),
                    ),
                    row=1,
                    col=i + 1,
                )
            fig.update_xaxes(title_text=xaxis + " [mmol/gdw/h]", row=1, col=i + 1)
            fig.update_yaxes(title_text=rxn[0] + " [mmol/gdw/h]", row=1, col=i + 1)

            fig.show()
