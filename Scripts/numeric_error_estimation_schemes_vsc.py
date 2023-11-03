
"""
Function library for the calculation of numeric control coefficients
- FCC method from Lu et al (2019)
- 1st order central finite difference calculation, approximate error: in the order of difference^2 (default 1e-4)
"""

def fcc_numeric_vsc_calculation(model, rxn, constraint_id, ref_obj, frac: float= 0.001):
    # change stoichiometry of UB constraint
    ref_stoich = model.constraints[constraint_id].get_linear_coefficients([rxn.forward_variable, rxn.reverse_variable])

    ref_fwd = ref_stoich[rxn.forward_variable]
    ref_rev = ref_stoich[rxn.reverse_variable]
    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: (1+frac) * ref_fwd,
        rxn.reverse_variable: (1+frac) * ref_rev
    })
    if ref_fwd==0 and ref_rev==0:
        return model, 0

    model.optimize()
    # if the model is infeasible, objective is assumed to be 0
    if model.solver.status == 'optimal':
        obj = model.objective.value
    else:
        obj = 0


    if ref_fwd !=0:
        FCC = (obj - ref_obj) / (ref_obj) / (((1+frac) * ref_fwd - ref_fwd) / ref_fwd)
    elif ref_rev !=0:
        FCC = -(obj - ref_obj) / (ref_obj) / (((1+frac) * ref_rev - ref_rev) / ref_rev)

    #reset the model
    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: ref_fwd,
        rxn.reverse_variable: ref_rev
    })
    return model, FCC

def fcc_numeric_vsc_optimizations(model, frac: float= 0.001):
    model.optimize()
    obj_value = model.objective.value  # v_z
    Crxn = []
    Cec = []
    past_enzymes = []
    #loop over the reactions and change the stoichiometry with 0.1% => optimize => calculate flux control coefficient
    for rxn in model.reactions:
        #change stoichiometry of UB constraint
        model, fcc_ub = fcc_numeric_vsc_calculation(model, rxn, f'{rxn.id}_ub', obj_value, frac)
        model, fcc_lb = fcc_numeric_vsc_calculation(model, rxn, f'{rxn.id}_lb', obj_value, frac)
        Crxn += [fcc_ub-fcc_lb]
        # check if reaction is associated with an enzyme:

        if not f'CE_{rxn.id}' in model.catalytic_events:
            continue
        enzymes = model.get_enzymes_with_reaction_id(rxn.id)
        for enzyme in enzymes:
            # perform analysis 1 time per enzyme
            if enzyme.id in past_enzymes:
                continue
            enz_var = model.enzyme_variables.get_by_id(enzyme.id)
            model, fcc_f = fcc_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_f', obj_value, frac)
            model, fcc_b = fcc_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_b', obj_value, frac)
            Cec += [fcc_f-fcc_b]
            past_enzymes.append(enzyme.id)

    return Crxn+ Cec

def fcc_numeric_evsc_optimizations(model, frac: float= 0.001):
    model.optimize()
    obj_value = model.objective.value  # v_z
    Cec = []
    #loop over the reactions and change the stoichiometry with 0.1% => optimize => calculate flux control coefficient
    for enzyme in model.enzymes:
        enz_var = model.enzyme_variables.get_by_id(enzyme.id)
        model, fcc_f = fcc_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_f', obj_value, frac)
        model, fcc_b = fcc_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_b', obj_value, frac)
        Cec += [fcc_f+fcc_b]
    return Cec


def first_forward_numeric_vsc_calculation(model, rxn, constraint_id, ref_obj, frac: float= 0.001):
    # change stoichiometry of UB constraint
    ref_stoich = model.constraints[constraint_id].get_linear_coefficients([rxn.forward_variable, rxn.reverse_variable])

    ref_fwd = ref_stoich[rxn.forward_variable]
    ref_rev = ref_stoich[rxn.reverse_variable]

    if ref_fwd==0 and ref_rev==0:
        return model, 0

    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: (1+frac) * ref_fwd,
        rxn.reverse_variable: (1+frac) * ref_rev
    })

    model.optimize()
    #if the model is infeasible, objective is assumed to be 0
    if model.solver.status == 'optimal':
        obj = model.objective.value
    else:
        obj = 0

    if ref_fwd !=0:
        FCC = (ref_obj - obj) / ((frac) * ref_fwd)
    elif ref_rev !=0:
        FCC = -(ref_obj - obj) / ((frac) * ref_rev)

    #reset the model
    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: ref_fwd,
        rxn.reverse_variable: ref_rev
    })
    return model, FCC

def first_forward_numeric_vsc_optimizations(model, frac: float= 0.001):
    model.optimize()
    obj_value = model.objective.value  # v_z
    Crxn = []
    Cec = []
    past_enzymes = []
    #loop over the reactions and change the stoichiometry with 0.1% => optimize => calculate flux control coefficient
    for rxn in model.reactions:
        #change stoichiometry of UB constraint
        model, fcc_ub = first_forward_numeric_vsc_calculation(model, rxn, f'{rxn.id}_ub', obj_value, frac)
        model, fcc_lb = first_forward_numeric_vsc_calculation(model, rxn, f'{rxn.id}_lb', obj_value, frac)
        Crxn += [fcc_ub-fcc_lb]

        #check if reaction is associated with an enzyme:
        if not f'CE_{rxn.id}' in model.catalytic_events:
            continue
        enzymes = model.get_enzymes_with_reaction_id(rxn.id)
        for enzyme in enzymes:
            # perform analysis 1 time per enzyme
            if enzyme.id in past_enzymes:
                continue
            enz_var = model.enzyme_variables.get_by_id(enzyme.id)
            model, fcc_f = first_forward_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_f', obj_value, frac)
            model, fcc_b = first_forward_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_b', obj_value, frac)
            Cec += [fcc_f-fcc_b]
            past_enzymes.append(enzyme.id)
    return Crxn+ Cec

def first_central_numeric_vsc_calculation(model, rxn, constraint_id, ref_obj, frac: float= 0.001):
    # change stoichiometry of UB constraint
    ref_stoich = model.constraints[constraint_id].get_linear_coefficients([rxn.forward_variable, rxn.reverse_variable])

    ref_fwd = ref_stoich[rxn.forward_variable]
    ref_rev = ref_stoich[rxn.reverse_variable]

    #if the variable is not associated, skip analysis
    if ref_fwd==0 and ref_rev==0:
        return model, 0

    if ref_obj==0:
        return model, 0

    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: (1+frac) * ref_fwd,
        rxn.reverse_variable: (1+frac) * ref_rev
    })
    model.optimize()
    # if the model is infeasible, objective is assumed to be 0
    if model.solver.status == 'optimal':
        obj_plus = model.objective.value
    else:
        obj_plus = 0

    ref_fwd = ref_stoich[rxn.forward_variable]
    ref_rev = ref_stoich[rxn.reverse_variable]
    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: (1 - frac) * ref_fwd,
        rxn.reverse_variable: (1 - frac) * ref_rev
    })
    model.optimize()
    # if the model is infeasible, objective is assumed to be 0
    if model.solver.status == 'optimal':
        obj_minus = model.objective.value
    else:
        obj_minus = 0

    if ref_fwd !=0:
        FCC = (obj_plus - obj_minus) / (2*(frac) * ref_fwd)  * ref_fwd/(ref_obj)
    elif ref_rev !=0:
        FCC = (obj_plus - obj_minus) / (2*(frac) * ref_rev) * ref_rev/(ref_obj)

    #reset the model
    model.constraints[constraint_id].set_linear_coefficients({
        rxn.forward_variable: ref_fwd,
        rxn.reverse_variable: ref_rev
    })
    return model, FCC

def first_central_numeric_vsc_optimizations(model, frac: float= 0.001):
    obj_value = model.objective.value  # v_z
    Crxn = []
    Cec = []
    past_enzymes = []
    #loop over the reactions and change the stoichiometry with 0.1% => optimize => calculate flux control coefficient
    for rxn in model.reactions:
        #change stoichiometry of UB constraint
        model, fcc_ub = first_central_numeric_vsc_calculation(model, rxn, f'{rxn.id}_ub', obj_value, frac)
        model, fcc_lb = first_central_numeric_vsc_calculation(model, rxn, f'{rxn.id}_lb', obj_value, frac)
        Crxn += [fcc_ub-fcc_lb]

        # check if reaction is associated with an enzyme:
        if not f'CE_{rxn.id}' in model.catalytic_events:
            continue
        enzymes = model.get_enzymes_with_reaction_id(rxn.id)
        for enzyme in enzymes:
            # perform analysis 1 time per enzyme
            if enzyme.id in past_enzymes:
                continue
            enz_var = model.enzyme_variables.get_by_id(enzyme.id)
        model, fcc_f = first_central_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_f', obj_value, frac)
        model, fcc_b = first_central_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_b', obj_value, frac)
        Cec += [fcc_f-fcc_b]
        past_enzymes.append(enzyme.id)

    return Crxn+Cec

def first_central_numeric_evsc_optimizations(model, frac: float= 0.001):
    obj_value = model.objective.value  # v_z
    Cec = []
    #loop over the reactions and change the stoichiometry with 0.1% => optimize => calculate flux control coefficient
    for enzyme in model.enzymes:
        enz_var = model.enzyme_variables.get_by_id(enzyme.id)
        model, fcc_f = first_central_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_f', obj_value, frac)
        model, fcc_b = first_central_numeric_vsc_calculation(model, enz_var, f'EC_{enzyme.id}_b', obj_value, frac)
        Cec += [fcc_f+fcc_b]

    return Cec