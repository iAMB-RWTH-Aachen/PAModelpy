import numpy as np
import pandas as pd
from warnings import warn
from typing import List, Optional, Union, Tuple, Dict

import logging
from optlang.symbolics import Zero
from optlang import Variable, interface

from cobra.core import Configuration, Model, Reaction
from cobra.util import solver as sutil
from cobra.util import ProcessPool
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.loopless import loopless_fva_iter

from ..Enzyme import Enzyme, EnzymeVariable


logger = logging.getLogger(__name__)
configuration = Configuration()

def _init_worker(model: "Model", loopless: bool, sense: str) -> None:
    """Initialize a global model object for multiprocessing.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    loopless: bool
        Whether to use loopless version.
    sense: {"max", "min"}
        Whether to maximise or minimise objective.

    """
    global _model
    global _loopless
    _model = model
    _model.solver.objective.direction = sense
    _loopless = loopless

def _fva_step(variable: Union[Variable, EnzymeVariable, Reaction]) -> Tuple[str, float]:
    """Take a step for calculating FVA.

    Parameters
    ----------
    variables: str
        The variable to test.

    Returns
    -------
    tuple of (str, float)
        The variable ID with the flux value.

    """
    global _model
    global _loopless
    # The previous objective assignment already triggers a reset
    # so directly update coefs here to not trigger redundant resets
    # in the history manager which can take longer than the actual
    # FVA for small models
    variable_id = variable.id if not (isinstance(variable, Variable) or not isinstance(variable, str)) else variable
    if isinstance(variable, Variable):
        model_var = _model.variables[variable.name]
        coeff = {model_var:1}
    elif variable_id in _model.variables:
        coeff = {_model.variables[variable_id]:1}
    else:
        model_var_f = _model.variables[variable.forward_variable.name]
        model_var_r = _model.variables[variable.reverse_variable.name]
        coeff = {model_var_f:1, model_var_r:-1}
    #set the fva objective
    _model.solver.objective.set_linear_coefficients(coeff)

    _model.slim_optimize()

    sutil.check_solver_status(_model.solver.status)
    if _loopless and isinstance(variable, Reaction):
        value = loopless_fva_iter(_model, variable)
    elif _loopless:
        warn('variable is not a reaction, thus cannot perform loopless FVA')
    else:
        value = _model.solver.objective.value

    # handle infeasible case
    if value is None:
        value = float("nan")
        logger.warning(
            f"Could not get a feasible solution for variable {variable}, setting it to NaN. "
            "This is usually due to numerical instability."
        )
    zero_coeff = {key:0 for key in coeff.keys()}
    _model.solver.objective.set_linear_coefficients(
       zero_coeff
    )
    return variable, value

def _get_variables(model, variable_type, variable_list, variable_mapping):
    """Collect (and optionally filter) model variables of a given *type*.

        This internal helper is used by the optimisation‑related utilities of the
        model.  Depending on the arguments it either

        1. **Returns all variables** of the requested ``variable_type`` (when
           ``variable_list`` is ``None``), **or**
        2. **Returns a subset** defined by ``variable_list`` – a mixture of
           variable IDs (``str``) and variable objects.

        The function also deals with *paired* forward/reverse variables that are
        stored as attributes ``forward_variable`` and ``reverse_variable`` on the
        supplied object (e.g. a reaction or an enzyme).

        Args:
            model: The modelling object that holds a ``variables`` dictionary and
                type‑specific containers (e.g. ``model.reactions``, ``model.enzymes``).
            variable_type (str or None): Key used to look up the correct attribute
                in ``variable_mapping`` (e.g. ``'reaction'`` → ``'reactions'``).
                If ``None`` the function will raise a ``KeyError`` when trying to
                access ``variable_mapping``.
            variable_list (list or None): List that defines which variables to
                return.  Elements may be:
                    * ``str`` – the exact variable ID,
                    * a COBRA‑style ``Variable`` object,
                    * an object that possesses ``id`` and optionally
                      ``forward_variable`` / ``reverse_variable`` attributes
                      (common for custom reaction/enzyme wrappers).
                If ``None`` all variables of the chosen ``variable_type`` are
                returned.
            variable_mapping (dict): Mapping from ``variable_type`` strings to the
                attribute name on ``model`` that stores the corresponding objects
                (e.g. ``{'reaction': 'reactions', 'enzyme': 'enzymes'}``).

        Returns:
            dict: Mapping ``variable_id → Variable`` containing the requested
            variables.  The dictionary is built by repeatedly merging sub‑dictionaries
            (using ``{**a, **b}``) so later entries overwrite earlier ones if IDs clash.

        Raises:
            KeyError: If ``variable_type`` is not present in ``variable_mapping``.
            AttributeError: If a non‑string entry in ``variable_list`` does not have
            the expected attributes (``id``, ``forward_variable``,
            ``reverse_variable``).

        Notes:
            * ``model.variables`` is assumed to behave like a dictionary that maps a
              variable’s ID (or name) to the actual variable object.
            * For reactions, the forward and reverse variables are retrieved
            * The function does **not** mutate the model; it only reads from it.
        """
    variables = {}
    if variable_list is not None:
        if variable_type is not None:
            variables = {var.id: var for var in getattr(model, variable_mapping[variable_type])}
        for var_id in variable_list:
            if isinstance(var_id, str):
                variables = {**variables, **{var_id: model.variables[var_id]}}
            else:
                if not isinstance(var_id, Variable):
                    if var_id.id in model.variables:
                        mapping = {var_id.id:model.variables[var_id.id]}
                    elif hasattr(var_id,'forward_variable'): #reactions and enzymes are associated with separate forward and reverse variables
                        mapping = {var.name: model.variables[var.name] for var in [var_id.forward_variable, var_id.reverse_variable]}
                else:
                    mapping = {var_id.name: model.variables[var_id.name]}
                variables = {**variables, **mapping}

    else:
        variables = {var.id: var for var in getattr(model, variable_mapping[variable_type])}

    return variables

def flux_variability_analysis(
    model: Model,
    variable_type: Optional[Union[Enzyme, EnzymeVariable, Reaction]] = None,
    variable_list: Optional[List[Union[Enzyme, EnzymeVariable, Reaction, str]]] = None,
    loopless: bool = False,
    fraction_of_optimum: float = 1.0,
    pfba_factor: Optional[float] = None,
    processes: Optional[int] = None,
        variable2attribute:Dict[Variable, str] = {
        Enzyme: 'enzyme_variables',
        EnzymeVariable: 'enzyme_variables',
        Reaction: 'reactions',
        }
) -> pd.DataFrame:
    """Determine the minimum and maximum flux value for each reaction.

    Adjusted from cobra.flux_analysis.variability.flux_variability_analysis

    Args:
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    variable_type: type of variable to perform FVA on
        if None, it will use all model variables
    variable_list : list of specified variables defined in variably_type or str, optional
        for which to obtain min/max fluxes. If None will use
        all of the specified variables in the model (default None). Will be added to all the variables
        defined in variable_type, or serve as total list of variables when variable_type is None
    loopless : bool, optional
        Whether to return only loopless solutions. This is significantly
        slower. Please also refer to the notes (default False).
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum (default 1.0).
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the `pfba_factor`
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds
        (default None).
    processes : int, optional
        The number of parallel processes to run. If not explicitly passed,
        will be set from the global configuration singleton (default None).
    variable2attribute : dict, optional
        Mapping between the variable types and the name of the attribute in
        which they are stored in the model.

    Returns
    pandas.DataFrame
        A data frame with variable identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    Notes:
        - This implements the fast version as described in [1]_. Please note that
        the flux distribution containing all minimal/maximal fluxes does not have
        to be a feasible solution for the model. Fluxes are minimized/maximized
        individually and a single minimal flux might require all others to be
        sub-optimal.

        - Using the loopless option will lead to a significant increase in
        computation time (about a factor of 100 for large models). However, the
        algorithm used here (see [2]_) is still more than 1000x faster than the
        'naive' version using `add_loopless(model)`. Also note that if you have
        included constraints that force a loop (for instance by setting all fluxes
        in a loop to be non-zero) this loop will be included in the solution.

    References:
    .. [1] Computationally efficient flux variability analysis.
       Gudmundsson S, Thiele I.
       BMC Bioinformatics. 2010 Sep 29;11:489.
       doi: 10.1186/1471-2105-11-489, PMID: 20920235

    .. [2] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions.
       Desouki AA, Jarre F, Gelius-Dietrich G, Lercher MJ.
       Bioinformatics. 2015 Jul 1;31(13):2159-65.
       doi: 10.1093/bioinformatics/btv096.

    """
    variables = _get_variables(model, variable_type, variable_list, variable2attribute)
    if processes is None:
        processes = configuration.processes

    num_variables = len(variables)
    processes = min(processes, num_variables)

    fva_result = pd.DataFrame(
        {
            "minimum": np.zeros(num_variables, dtype=float),
            "maximum": np.zeros(num_variables, dtype=float),
        },
        index=variables,
    )
    prob = model.problem
    with model:
        # Safety check before setting up FVA.
        model.slim_optimize(
            error_value=None,
            message="There is no optimal solution for the chosen objective!",
        )
        # Add the previous objective as a variable to the model then set it to
        # zero. This also uses the fraction to create the lower/upper bound for
        # the old objective.
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value,
            )
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value,
            )
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.0:
                warn(
                    "The 'pfba_factor' should be larger or equal to 1.",
                    UserWarning,
                )
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum,
                    lb=0,
                    ub=0,
                    name="flux_sum_constraint",
                )
            model.add_cons_vars([flux_sum, flux_sum_constraint])

        model.objective = Zero  # This will trigger the reset as well
        for what in ("minimum", "maximum"):
            if processes > 1:
                # We create and destroy a new pool here in order to set the
                # objective direction for all reactions. This creates a
                # slight overhead but seems the most clean.
                chunk_size = len(variables) // processes
                with ProcessPool(
                    processes,
                    initializer=_init_worker,
                    initargs=(model, loopless, what[:3]),
                ) as pool:
                    for var, value in pool.imap_unordered(
                        _fva_step, variables.values(), chunksize=chunk_size
                    ):
                        try: var_id = var.id
                        except: var_id = var.name
                        fva_result.at[var_id, what] = value
            else:
                _init_worker(model, loopless, what[:3])
                for var, value in map(_fva_step, variables):
                    try: var_id = var.id
                    except: var_id = var.name
                    fva_result.at[var_id, what] = value

    return fva_result[["minimum", "maximum"]]