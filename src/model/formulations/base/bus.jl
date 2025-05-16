# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

function _add_bus!(
    model::JuMP.Model,
    b::Bus,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    net_injection = _init(model, :expr_net_injection)
    curtail = _init(model, :curtail)
    bus_load = model[:bus_load]
    for t in 1:nCont+nInt
        # Fixed load
        # net_injection[sc.name, b.name, t] = AffExpr(-b.load[t])
        bus_load[b.name, t] = @variable(model, base_name="bus_load_$(b.name)_$(t)")
        net_injection[sc.name, b.name, t] = AffExpr()
        add_to_expression!(
            net_injection[sc.name, b.name, t],
            bus_load[b.name, t],
            -1.0,
        )

        # Load curtailment
        curtail[sc.name, b.name, t] =
            @variable(model, lower_bound = 0, upper_bound = b.load[t], base_name = "curtail_$(b.name)_$(t)")

        add_to_expression!(
            net_injection[sc.name, b.name, t],
            curtail[sc.name, b.name, t],
            1.0,
        )
        add_to_expression!(
            model[:obj],
            curtail[sc.name, b.name, t],
            sc.power_balance_penalty[t] * sc.probability,
        )
    end
    return
end
