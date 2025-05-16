# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

function _add_status_vars!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_status_vars::Gar1962.StatusVars, 
    nCont::Int, 
    nInt::Int 
)::Nothing
    is_on = _init(model, :is_on)
    switch_on = _init(model, :switch_on)
    switch_off = _init(model, :switch_off)
    initial_on = _init(model, :initial_on) # if g.initial_status > 0, initial_on = 1.0
    initial_status = _init(model, :initial_status)
    initial_on[g.name] = @variable(model, binary = true, base_name = "initial_on_$(g.name)")
    initial_status[g.name] = @variable(model, binary = false, base_name = "initial_status_$(g.name)")
    for t in 1:nCont+nInt
        # if g.must_run[t]
        #     printstyled("$(g.name) must be on at time $(t)\n"; color=:red)
        #     is_on[g.name, t] = 1.0
        #     # switch_on[g.name, t] = (t == 1 ? 1.0 - _is_initially_on(g) : 0.0)
        #     switch_on[g.name, t] = (t == 1 ? 1.0 - initial_on[g.name] : 0.0)
        #     switch_off[g.name, t] = 0.0
        if t <= nInt
            is_on[g.name, t] = @variable(model, binary = true, base_name = "ison_int_$(g.name)_$(t)")
            switch_on[g.name, t] = @variable(model, binary = true, base_name = "switchon_int_$(g.name)_$(t)")
            switch_off[g.name, t] = @variable(model, binary = true, base_name = "switchoff_int_$(g.name)_$(t)")
        else
            is_on[g.name, t] = @variable(model, binary = false,lower_bound=0.0,upper_bound=1.0, base_name = "ison_cont_$(g.name)_$(t)")
            switch_on[g.name, t] = @variable(model, binary = false,lower_bound=0.0,upper_bound=1.0, base_name = "switchon_cont_$(g.name)_$(t)")
            switch_off[g.name, t] = @variable(model, binary = false,lower_bound=0.0,upper_bound=1.0, base_name = "switchoff_cont_$(g.name)_$(t)")
        end
        add_to_expression!(model[:obj], is_on[g.name, t], g.min_power_cost[t]) 
    end
    return
end


function _update_status_vars!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_status_vars::Gar1962.StatusVars, 
    nCont::Int, 
    nInt::Int,
    offset::Int
)::Nothing
    is_on = _init(model, :is_on)
    switch_on = _init(model, :switch_on)
    switch_off = _init(model, :switch_off)
    for t in 1+offset:offset+nCont+nInt
        if g.must_run[t]
            is_on[g.name, t-offset] = 1.0
            switch_on[g.name, t-offset] = (t == 1 ? 1.0 - _is_initially_on(g) : 0.0)
            switch_off[g.name, t-offset] = 0.0
        end
        # No need to add vars
        # elseif t <= nInt
        #     is_on[g.name, t] = @variable(model, binary = true, base_name = "ison_int_$(g.name)_$(t)")
        #     switch_on[g.name, t] = @variable(model, binary = true, base_name = "switchon_int_$(g.name)_$(t)")
        #     switch_off[g.name, t] = @variable(model, binary = true, base_name = "switchoff_int_$(g.name)_$(t)")
        # else
        #     is_on[g.name, t] = @variable(model, binary = false,lower_bound=0.0,upper_bound=1.0, base_name = "ison_cont_$(g.name)_$(t)")
        #     switch_on[g.name, t] = @variable(model, binary = false,lower_bound=0.0,upper_bound=1.0, base_name = "switchon_cont_$(g.name)_$(t)")
        #     switch_off[g.name, t] = @variable(model, binary = false,lower_bound=0.0,upper_bound=1.0, base_name = "switchoff_cont_$(g.name)_$(t)")
        # end
        add_to_expression!(model[:obj], is_on[g.name, t-offset], g.min_power_cost[t]) 
    end
    return
end

function _add_status_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_status_vars::Gar1962.StatusVars, 
    nCont::Int, 
    nInt::Int
)::Nothing
    eq_binary_link = _init(model, :eq_binary_link)
    eq_switch_on_off = _init(model, :eq_switch_on_off)
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]
    initial_on = model[:initial_on]
    for t in 1:nCont+nInt
        if !g.must_run[t]
            # Link binary variables
            if t == 1
                # eq_binary_link[g.name, t] = @constraint(
                #     model,
                #     is_on[g.name, t] - _is_initially_on(g) ==
                #     switch_on[g.name, t] - switch_off[g.name, t],
                #     base_name = "eq_binary_link_$(g.name)_$(t)"
                # )
                eq_binary_link[g.name, t] = @constraint(
                    model,
                    is_on[g.name, t] - initial_on[g.name] ==
                    switch_on[g.name, t] - switch_off[g.name, t],
                    base_name = "eq_binary_link_$(g.name)_$(t)"
                )
            else
                eq_binary_link[g.name, t] = @constraint(
                    model,
                    is_on[g.name, t] - is_on[g.name, t-1] ==
                    switch_on[g.name, t] - switch_off[g.name, t],
                    base_name = "eq_binary_link_$(g.name)_$(t)"
                )
            end
            # Cannot switch on and off at the same time
            eq_switch_on_off[g.name, t] = @constraint(
                model,
                switch_on[g.name, t] + switch_off[g.name, t] <= 1,
                base_name = "eq_switch_on_off_$(g.name)_$(t)"
            )
        end
    end
    return
end
