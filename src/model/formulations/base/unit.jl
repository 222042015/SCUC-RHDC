# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

# Function for adding variables, constraints, and objective function terms
# related to the binary commitment, startup and shutdown decisions of units
function _add_unit_commitment!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation::Formulation,
    nCont::Int,
    nInt::Int
)
    if !all(g.must_run) && any(g.must_run)
        error("Partially must-run units are not currently supported")
    end
    if g.initial_power === nothing || g.initial_status === nothing
        error("Initial conditions for $(g.name) must be provided")
    end

    # Variables
    _add_startup_shutdown_vars!(model, g, nCont, nInt) #DONE
    _add_status_vars!(model, g, formulation.status_vars, nCont, nInt) #DONE

    # Constraints and objective function
    _add_min_uptime_downtime_eqs!(model, g, nCont, nInt) #DONE
    _add_startup_cost_eqs!(model, g, formulation.startup_costs, nCont, nInt) #TODO:: startup costs
    _add_status_eqs!(model, g, formulation.status_vars, nCont, nInt) # DONE
    _add_commitment_status_eqs!(model, g, nCont, nInt) # DONE::not involved
    return
end


# Function for adding variables, constraints, and objective function terms
# related to the continuous dispatch decisions of units
function _add_unit_dispatch!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation::Formulation,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)

    # Variables
    _add_production_vars!(model, g, formulation.prod_vars, sc, nCont, nInt)
    _add_spinning_reserve_vars!(model, g, sc, nCont, nInt)
    _add_flexiramp_reserve_vars!(model, g, sc, nCont, nInt)

    # Constraints and objective function
    _add_net_injection_eqs!(model, g, sc, nCont, nInt)
    _add_production_limit_eqs!(model, g, formulation.prod_vars, sc, nCont, nInt) # DONE
    _add_production_piecewise_linear_eqs!( # DONE
        model,
        g,
        formulation.prod_vars,
        formulation.pwl_costs,
        formulation.status_vars,
        sc,
        nCont, 
        nInt
    )
    _add_ramp_eqs!( # DONE
        model,
        g,
        formulation.prod_vars,
        formulation.ramping,
        formulation.status_vars,
        sc,
        nCont, 
        nInt
    )
    _add_startup_shutdown_limit_eqs!(model, g, sc, nCont, nInt) #DONE
    return
end

# _is_initially_on(g::ThermalUnit)::Float64 = (g.initial_status > 0 ? 1.0 : 0.0)

function _add_spinning_reserve_vars!(
    model::JuMP.Model,
    g::ThermalUnit,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    # same, all cont var, no need to change according to time
    reserve = _init(model, :reserve)
    total_reserve = _init(model, :total_reserve)
    reserve_shortfall = _init(model, :reserve_shortfall)
    for r in g.reserves
        r.type == "spinning" || continue
        total_reserve[sc.name, g.name, 0] = @variable(model, base_name="total_reserve_$(g.name)_0")
        for t in 1:nCont+nInt
            total_reserve[sc.name, g.name, t] = @variable(model, base_name="total_reserve_$(g.name)_$(t)")
            reserve[sc.name, r.name, g.name, t] =
                @variable(model, lower_bound = 0, base_name="reserve_$(g.name)_$(t)")
            if (sc.name, r.name, t) ∉ keys(reserve_shortfall)
                reserve_shortfall[sc.name, r.name, t] =
                    @variable(model, lower_bound = 0, base_name="reserve_shortfall_$(g.name)_$(t)")
                if r.shortfall_penalty < 0
                    set_upper_bound(reserve_shortfall[sc.name, r.name, t], 0.0)
                end
            end
        end
    end
    return
end

function _add_flexiramp_reserve_vars!(
    model::JuMP.Model,
    g::ThermalUnit,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    upflexiramp = _init(model, :upflexiramp)
    upflexiramp_shortfall = _init(model, :upflexiramp_shortfall)
    mfg = _init(model, :mfg)
    dwflexiramp = _init(model, :dwflexiramp)
    dwflexiramp_shortfall = _init(model, :dwflexiramp_shortfall)
    for t in 1:nCont+nInt
        # maximum feasible generation, \bar{g_{its}} in Wang & Hobbs (2016)
        mfg[sc.name, g.name, t] = @variable(model, lower_bound = 0, base_name = "mfg_$(g.name)_$(t)")
        for r in g.reserves
            r.type == "flexiramp" || continue
            upflexiramp[sc.name, r.name, g.name, t] = @variable(model) # up-flexiramp, ur_{it} in Wang & Hobbs (2016)
            dwflexiramp[sc.name, r.name, g.name, t] = @variable(model) # down-flexiramp, dr_{it} in Wang & Hobbs (2016)
            if (sc.name, r.name, t) ∉ keys(upflexiramp_shortfall)
                upflexiramp_shortfall[sc.name, r.name, t] =
                    @variable(model, lower_bound = 0, base_name="upflexramp_$(r.name)_$(t)")
                dwflexiramp_shortfall[sc.name, r.name, t] =
                    @variable(model, lower_bound = 0, base_name="dwflexramp_$(r.name)_$(t)")
                if r.shortfall_penalty < 0
                    set_upper_bound(
                        upflexiramp_shortfall[sc.name, r.name, t],
                        0.0,
                    )
                    set_upper_bound(
                        dwflexiramp_shortfall[sc.name, r.name, t],
                        0.0,
                    )
                end
            end
        end
    end
    return
end

function _add_startup_shutdown_vars!(model::JuMP.Model, g::ThermalUnit, nCont::Int, nInt::Int )::Nothing
    startup = _init(model, :startup)
    switch_off_fixed = _init(model, :switch_off_fixed)
    initial_sum = _init(model, :initial_sum)
    for s in 1:length(g.startup_categories)
        for t in 1:nInt
            startup[g.name, t, s] = @variable(model, binary = true, base_name = "startup_int_$(g.name)_$(t)" )
            if s == length(g.startup_categories)
                continue
            end
            initial_sum[g.name, t, s] = @variable(model, binary = true, base_name = "initial_sum_$(g.name)_$(t)_$(s)" )
        end
        for t in nInt+1:nInt+nCont
            startup[g.name, t, s] = @variable(model, binary = false, lower_bound=0.0,upper_bound=1.0, base_name = "startup_cont_$(g.name)_$(t)" )
            if s == length(g.startup_categories)
                continue
            end
            initial_sum[g.name, t, s] = @variable(model, binary = true, base_name = "initial_sum_$(g.name)_$(t)_$(s)" )
        end
    end


    for t in 1:nCont+nInt
        for s in 1:length(g.startup_categories)-1
            range_start = t - g.startup_categories[s+1].delay + 1
            range_end = t - g.startup_categories[s].delay
            range = (range_start:range_end)
            for i in range
                if i < 1
                    switch_off_fixed[g.name, i, s] = @variable(model, binary = true, base_name = "switch_off_fixed_$(g.name)_$(t)_$(s)" )
                end
            end
        end
    end

    return
end

function _add_startup_shutdown_limit_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    eq_shutdown_limit = _init(model, :eq_shutdown_limit)
    eq_startup_limit = _init(model, :eq_startup_limit)
    is_on = model[:is_on]
    prod_above = model[:prod_above]
    # reserve = _total_reserves(model, g, sc, nCont, nInt)
    total_reserve = model[:total_reserve]
    # spinning_reserves = [r for r in g.reserves if r.type == "spinning"]
    # if !isempty(spinning_reserves)
    #     for t in 1:nCont+nInt
    #         total_reserve[sc.name, g.name, t] = sum(model[:reserve][sc.name, r.name, g.name, t] for r in spinning_reserves)
    #     end
    # else
    #     for t in 1:nCont+nInt
    #         total_reserve[sc.name, g.name, t] = 0.0
    #     end
    # end
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]
    T = nCont+nInt
    for t in 1:T
        # Startup limit
        eq_startup_limit[sc.name, g.name, t] = @constraint(
            model,
            prod_above[sc.name, g.name, t] + total_reserve[sc.name, g.name, t] <=
            (g.max_power[t] - g.min_power[t]) * is_on[g.name, t] -
            max(0, g.max_power[t] - g.startup_limit) * switch_on[g.name, t],
            base_name = "eq_startup_limit_$(sc.name)_$(g.name)_$(t)"
        )
        # # Shutdown limit
        # if t == 1 && g.initial_power > g.shutdown_limit
        #     eq_shutdown_limit[sc.name, g.name, 0] =
        #         @constraint(model, switch_off[g.name, 1] <= 0, 
        #             base_name="eq_shutdown_limit_$(sc.name)_$(g.name)_0")
        # end
        if t < T
            eq_shutdown_limit[sc.name, g.name, t] = @constraint(
                model,
                prod_above[sc.name, g.name, t] <=
                (g.max_power[t] - g.min_power[t]) * is_on[g.name, t] -
                max(0, g.max_power[t] - g.shutdown_limit) *
                switch_off[g.name, t+1],
                base_name = "eq_shutdown_limit_$(sc.name)_$(g.name)_$(t)"
            )
        end
    end
    return
end

function _add_ramp_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation::RampingFormulation,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    prod_above = model[:prod_above]
    reserve = _total_reserves(model, g, sc, nCont, nInt)
    eq_ramp_up = _init(model, :eq_ramp_up)
    eq_ramp_down = _init(model, :eq_ramp_down)
    for t in 1:nCont+nInt
        # Ramp up limit 
        if t == 1
            if _is_initially_on(g) == 1
                eq_ramp_up[sc.name, g.name, t] = @constraint(
                    model,
                    prod_above[sc.name, g.name, t] + reserve[t] <=
                    (g.initial_power - g.min_power[t]) + g.ramp_up_limit
                )
            end
        else
            eq_ramp_up[sc.name, g.name, t] = @constraint(
                model,
                prod_above[sc.name, g.name, t] + reserve[t] <=
                prod_above[sc.name, g.name, t-1] + g.ramp_up_limit
                # , base_name="ramp_up_limit_$(g.name)_$(t)"
            )
        end

        # Ramp down limit
        if t == 1
            if _is_initially_on(g) == 1
                eq_ramp_down[sc.name, g.name, t] = @constraint(
                    model,
                    prod_above[sc.name, g.name, t] >=
                    (g.initial_power - g.min_power[t]) - g.ramp_down_limit
                )
            end
        else
            eq_ramp_down[sc.name, g.name, t] = @constraint(
                model,
                prod_above[sc.name, g.name, t] >=
                prod_above[sc.name, g.name, t-1] - g.ramp_down_limit
            )
        end
    end
end

function _add_min_uptime_downtime_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    nCont::Int,
    nInt::Int
)::Nothing
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]
    eq_min_uptime = _init(model, :eq_min_uptime)
    eq_min_downtime = _init(model, :eq_min_downtime)
    T = nCont+nInt
    for t in 1:nCont+nInt
        # Minimum up-time
        eq_min_uptime[g.name, t] = @constraint(
            model,
            sum(switch_on[g.name, i] for i in (t-g.min_uptime+1):t if i >= 1) <= is_on[g.name, t],
            base_name="eq_min_uptime_$(g.name)_$(t)"
        )
        # Minimum down-time
        eq_min_downtime[g.name, t] = @constraint(
            model,
            sum(
                switch_off[g.name, i] for i in (t-g.min_downtime+1):t if i >= 1
            ) <= 1 - is_on[g.name, t],
            base_name="eq_min_downtime_$(g.name)_$(t)"
        )
        # # Minimum up/down-time for initial periods
        # if t == 1
        #     if g.initial_status > 0
        #         eq_min_uptime[g.name, 0] = @constraint(
        #             model,
        #             sum(
        #                 switch_off[g.name, i] for
        #                 i in 1:(g.min_uptime-g.initial_status) if i <= T
        #             ) == 0,
        #             base_name="eq_min_uptime_$(g.name)_0"
        #         )
        #     else
        #         eq_min_downtime[g.name, 0] = @constraint(
        #             model,
        #             sum(
        #                 switch_on[g.name, i] for
        #                 i in 1:(g.min_downtime+g.initial_status) if i <= T
        #             ) == 0,
        #             base_name="eq_min_downtime_$(g.name)_0"
        #         )
        #     end
        # end
    end
end


# NEED PREVIOUS STARTUP SHUTDOWN STATUS
function _update_min_uptime_downtime_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    nCont::Int,
    nInt::Int,
    offset::Int
)::Nothing
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]
    eq_min_uptime = _init(model, :eq_min_uptime)
    eq_min_downtime = _init(model, :eq_min_downtime)
    T = nCont+nInt
    for t in 1:nCont+nInt
        # Minimum up-time
        eq_min_uptime[g.name, t] = @constraint(
            model,
            sum(switch_on[g.name, i] for i in (t-g.min_uptime+1):t if i >= 1) <= is_on[g.name, t]
        )
        # Minimum down-time
        eq_min_downtime[g.name, t] = @constraint(
            model,
            sum(
                switch_off[g.name, i] for i in (t-g.min_downtime+1):t if i >= 1
            ) <= 1 - is_on[g.name, t]
        )
        # # Minimum up/down-time for initial periods
        # if t == 1
        #     if g.initial_status > 0
        #         eq_min_uptime[g.name, 0] = @constraint(
        #             model,
        #             sum(
        #                 switch_off[g.name, i] for
        #                 i in 1:(g.min_uptime-g.initial_status) if i <= T
        #             ) == 0
        #         )
        #     else
        #         eq_min_downtime[g.name, 0] = @constraint(
        #             model,
        #             sum(
        #                 switch_on[g.name, i] for
        #                 i in 1:(g.min_downtime+g.initial_status) if i <= T
        #             ) == 0
        #         )
        #     end
        # end
    end
end

function _add_commitment_status_eqs!(model::JuMP.Model, g::ThermalUnit, nCont::Int, nInt::Int)::Nothing
    is_on = model[:is_on]
    T = nCont+nInt
    eq_commitment_status = _init(model, :eq_commitment_status)
    for t in 1:T
        if g.commitment_status[t] !== nothing
            eq_commitment_status[g.name, t] = @constraint(
                model,
                is_on[g.name, t] == (g.commitment_status[t] ? 1.0 : 0.0)
            )
        end
    end
    return
end

function _add_net_injection_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    expr_net_injection = model[:expr_net_injection]
    for t in 1:nCont+nInt
        # Add to net injection expression
        add_to_expression!(
            expr_net_injection[sc.name, g.bus.name, t],
            model[:prod_above][sc.name, g.name, t],
            1.0,
        )
        add_to_expression!(
            expr_net_injection[sc.name, g.bus.name, t],
            model[:is_on][g.name, t],
            g.min_power[t],
        )
    end
end

function _total_reserves(model, g, sc, nCont, nInt)::Vector
    T = nCont+nInt
    reserve = [0.0 for _ in 1:T]
    spinning_reserves = [r for r in g.reserves if r.type == "spinning"]
    if !isempty(spinning_reserves)
        reserve += [
            sum(
                model[:reserve][sc.name, r.name, g.name, t] for
                r in spinning_reserves
            ) for t in 1:T
        ]
    end
    return reserve
end