# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

using JuMP

"""
    generate_initial_conditions!(instance, optimizer)

Generates feasible initial conditions for the given instance, by constructing
and solving a single-period mixed-integer optimization problem, using the given
optimizer. The instance is modified in-place.
"""
function generate_initial_conditions!(
    instance::ucRHInstance,
    optimizer,
)::Nothing
    # Process first scenario
    _generate_initial_conditions!(instance.scenarios[1], optimizer)

    # Copy initial conditions to remaining scenarios
    for (si, sc) in enumerate(instance.scenarios)
        si > 1 || continue
        for (gi, g) in sc.thermal_units
            g_ref = instance.scenarios[1].thermal_units[gi]
            g.initial_power = g_ref.initial_power
            g.initial_status = g_ref.initial_status
        end
    end
end

function _generate_initial_conditions!(
    sc::ucRHScenario,
    optimizer,
)::Nothing
    G = sc.thermal_units
    B = sc.buses
    PU = sc.profiled_units
    t = 1
    mip = JuMP.Model(optimizer)

    # Decision variables
    @variable(mip, x[G], Bin)
    @variable(mip, p[G] >= 0)
    @variable(mip, pu[PU])

    # Constraint: Minimum power
    @constraint(mip, min_power[g in G], p[g] >= g.min_power[t] * x[g])
    @constraint(mip, pu_min_power[k in PU], pu[k] >= k.min_power[t])

    # Constraint: Maximum power
    @constraint(mip, max_power[g in G], p[g] <= g.max_power[t] * x[g])
    @constraint(mip, pu_max_power[k in PU], pu[k] <= k.max_power[t])

    # Constraint: Production equals demand
    @constraint(
        mip,
        power_balance,
        sum(b.load[t] for b in B) ==
        sum(p[g] for g in G) + sum(pu[k] for k in PU)
    )

    # Constraint: Must run
    for g in G
        if g.must_run[t]
            @constraint(mip, x[g] == 1)
        end
    end

    # Objective function
    function cost_slope(g)
        mw = g.min_power[t]
        c = g.min_power_cost[t]
        for k in g.cost_segments
            mw += k.mw[t]
            c += k.mw[t] * k.cost[t]
        end
        if mw < 1e-3
            return 0.0
        else
            return c / mw
        end
    end
    @objective(
        mip,
        Min,
        sum(p[g] * cost_slope(g) for g in G) +
        sum(pu[k] * k.cost[t] for k in PU)
    )

    JuMP.optimize!(mip)

    for g in G
        if JuMP.value(x[g]) > 0
            g.initial_power = JuMP.value(p[g])
            g.initial_status = 24
        else
            g.initial_power = 0
            g.initial_status = -24
        end
    end
    return
end
