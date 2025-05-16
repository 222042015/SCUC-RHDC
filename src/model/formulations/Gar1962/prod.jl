# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

function _add_production_vars!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_prod_vars::Gar1962.ProdVars,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    prod_above = _init(model, :prod_above)
    segprod = _init(model, :segprod)
    # all cont var, no need to change type]
    prod_above[sc.name, g.name, 0] = @variable(model, lower_bound = 0, base_name = "prodabv_$(g.name)_0")
    for t in 1:nCont+nInt
        for k in 1:length(g.cost_segments)
            segprod[sc.name, g.name, t, k] = @variable(model, lower_bound = 0, base_name = "segprod_$(g.name)_$(t)")
        end
        prod_above[sc.name, g.name, t] = @variable(model, lower_bound = 0, base_name = "prodabv_$(g.name)_$(t)")
    end
    return
end

function _add_production_limit_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_prod_vars::Gar1962.ProdVars,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    eq_prod_limit = _init(model, :eq_prod_limit)
    is_on = model[:is_on]
    prod_above = model[:prod_above]
    # reserve = _total_reserves(model, g, sc, nCont, nInt)
    total_reserve = model[:total_reserve]
    spinning_reserves = [r for r in g.reserves if r.type == "spinning"]
    if !isempty(spinning_reserves)
        for t in 1:nCont+nInt
            total_reserve[sc.name, g.name, t] = sum(model[:reserve][sc.name, r.name, g.name, t] for r in spinning_reserves)
        end
    else
        for t in 1:nCont+nInt
            total_reserve[sc.name, g.name, t] = 0.0
        end
    end
    gn = g.name
    for t in 1:nCont+nInt
        # Objective function terms for production costs
        # Part of (69) of Kneuven et al. (2020) as C^R_g * u_g(t) term

        # Production limit
        # Equation (18) in Kneuven et al. (2020)
        #   as \bar{p}_g(t) \le \bar{P}_g u_g(t)
        # amk: this is a weaker version of (20) and (21) in Kneuven et al. (2020)
        #      but keeping it here in case those are not present
        power_diff = max(g.max_power[t], 0.0) - max(g.min_power[t], 0.0)
        if power_diff < 1e-7
            power_diff = 0.0
        end
        eq_prod_limit[sc.name, gn, t] = @constraint(
            model,
            # prod_above[sc.name, gn, t] + reserve[t] <= power_diff * is_on[gn, t],
            prod_above[sc.name, gn, t] + total_reserve[sc.name, gn, t] <= power_diff * is_on[gn, t],
            base_name = "eq_prod_limit_$(gn)_$(t)"
        )
    end
end
