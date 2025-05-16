# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

# function _add_ramp_eqs!(
#     model::JuMP.Model,
#     g::ThermalUnit,
#     formulation_prod_vars::Gar1962.ProdVars,
#     formulation_ramping::MorLatRam2013.Ramping,
#     formulation_status_vars::Gar1962.StatusVars,
#     sc::ucRHScenario,
#     nCont::Int,
#     nInt::Int
# )::Nothing
#     # TODO: Move upper case constants to model[:instance]
#     RESERVES_WHEN_START_UP = true
#     RESERVES_WHEN_RAMP_UP = true
#     RESERVES_WHEN_RAMP_DOWN = true
#     RESERVES_WHEN_SHUT_DOWN = true
#     is_initially_on = (g.initial_status > 0)
#     SU = g.startup_limit
#     SD = g.shutdown_limit
#     RU = g.ramp_up_limit
#     RD = g.ramp_down_limit
#     gn = g.name
#     eq_ramp_down = _init(model, :eq_ramp_down)
#     eq_ramp_up = _init(model, :eq_str_ramp_up)
#     reserve = _total_reserves(model, g, sc, nCont, nInt)

#     # Gar1962.ProdVars
#     prod_above = model[:prod_above]

#     # Gar1962.StatusVars
#     is_on = model[:is_on]
#     switch_off = model[:switch_off]
#     switch_on = model[:switch_on]

#     for t in 1:nCont+nInt
#         time_invariant =
#             (t > 1) ? (abs(g.min_power[t] - g.min_power[t-1]) < 1e-7) : true

#         # Ramp up limit
#         if t == 1
#             if is_initially_on
#                 eq_ramp_up[sc.name, gn, t] = @constraint(
#                     model,
#                     g.min_power[t] +
#                     prod_above[sc.name, gn, t] +
#                     (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) <=
#                     g.initial_power + RU
#                 )
#             end
#         else
#             # amk: without accounting for time-varying min power terms,
#             #      we might get an infeasible schedule, e.g. if min_power[t-1] = 0, min_power[t] = 10
#             #      and ramp_up_limit = 5, the constraint (p'(t) + r(t) <= p'(t-1) + RU)
#             #      would be satisfied with p'(t) = r(t) = p'(t-1) = 0
#             #      Note that if switch_on[t] = 1, then eqns (20) or (21) go into effect
#             if !time_invariant
#                 # Use equation (24) instead
#                 SU = g.startup_limit
#                 max_prod_this_period =
#                     g.min_power[t] * is_on[gn, t] +
#                     prod_above[sc.name, gn, t] +
#                     (
#                         RESERVES_WHEN_START_UP || RESERVES_WHEN_RAMP_UP ?
#                         reserve[t] : 0.0
#                     )
#                 min_prod_last_period =
#                     g.min_power[t-1] * is_on[gn, t-1] +
#                     prod_above[sc.name, gn, t-1]
#                 eq_ramp_up[gn, t] = @constraint(
#                     model,
#                     max_prod_this_period - min_prod_last_period <=
#                     RU * is_on[gn, t-1] + SU * switch_on[gn, t]
#                 )
#             else
#                 # Equation (26) in Kneuven et al. (2020)
#                 # TODO: what if RU < SU? places too stringent upper bound
#                 # prod_above[gn, t] when starting up, and creates diff with (24).
#                 eq_ramp_up[sc.name, gn, t] = @constraint(
#                     model,
#                     prod_above[sc.name, gn, t] +
#                     (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) -
#                     prod_above[sc.name, gn, t-1] <= RU
#                 )
#             end
#         end

#         # Ramp down limit
#         if t == 1
#             if is_initially_on
#                 # TODO If RD < SD, or more specifically if
#                 #        min_power + RD < initial_power < SD
#                 #      then the generator should be able to shut down at time t = 1,
#                 #      but the constraint below will force the unit to produce power
#                 eq_ramp_down[sc.name, gn, t] = @constraint(
#                     model,
#                     g.initial_power -
#                     (g.min_power[t] + prod_above[sc.name, gn, t]) <= RD
#                 )
#             end
#         else
#             # amk: similar to ramp_up, need to account for time-dependent min_power
#             if !time_invariant
#                 # Revert to (25)
#                 SD = g.shutdown_limit
#                 max_prod_last_period =
#                     g.min_power[t-1] * is_on[gn, t-1] +
#                     prod_above[sc.name, gn, t-1] +
#                     (
#                         RESERVES_WHEN_SHUT_DOWN || RESERVES_WHEN_RAMP_DOWN ?
#                         reserve[t-1] : 0.0
#                     )
#                 min_prod_this_period =
#                     g.min_power[t] * is_on[gn, t] + prod_above[sc.name, gn, t]
#                 eq_ramp_down[gn, t] = @constraint(
#                     model,
#                     max_prod_last_period - min_prod_this_period <=
#                     RD * is_on[gn, t] + SD * switch_off[gn, t]
#                 )
#             else
#                 # Equation (27) in Kneuven et al. (2020)
#                 # TODO: Similar to above, what to do if shutting down in time t
#                 # and RD < SD? There is a difference with (25).
#                 eq_ramp_down[sc.name, gn, t] = @constraint(
#                     model,
#                     prod_above[sc.name, gn, t-1] +
#                     (RESERVES_WHEN_RAMP_DOWN ? reserve[t-1] : 0.0) -
#                     prod_above[sc.name, gn, t] <= RD
#                 )
#             end
#         end
#     end
# end


function _add_ramp_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_prod_vars::Gar1962.ProdVars,
    formulation_ramping::MorLatRam2013.Ramping,
    formulation_status_vars::Gar1962.StatusVars,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_START_UP = true
    RESERVES_WHEN_RAMP_UP = true
    RESERVES_WHEN_RAMP_DOWN = true
    RESERVES_WHEN_SHUT_DOWN = true
    is_initially_on = (g.initial_status > 0)
    SU = g.startup_limit
    SD = g.shutdown_limit
    RU = g.ramp_up_limit
    RD = g.ramp_down_limit
    gn = g.name
    eq_ramp_down = _init(model, :eq_ramp_down)
    eq_ramp_up = _init(model, :eq_str_ramp_up)
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


    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    for t in 1:nCont+nInt
        # time_invariant =
        #     (t > 1) ? (abs(g.min_power[t] - g.min_power[t-1]) < 1e-7) : true

        # eq_ramp_up[sc.name, gn, t] = @constraint(
        #         model,
        #         prod_above[sc.name, gn, t] +
        #         (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) -
        #         prod_above[sc.name, gn, t-1] <= RU,
        #         base_name = "eq_ramp_up_$(sc.name)_$(gn)_$(t)"
        # )

        # eq_ramp_down[sc.name, gn, t] = @constraint(
        #             model,
        #             prod_above[sc.name, gn, t-1] +
        #             (RESERVES_WHEN_RAMP_DOWN ? reserve[t-1] : 0.0) -
        #             prod_above[sc.name, gn, t] <= RD,
        #             base_name = "eq_ramp_dn_$(sc.name)_$(gn)_$(t)"
        #         )

        eq_ramp_up[sc.name, gn, t] = @constraint(
                model,
                prod_above[sc.name, gn, t] +
                (RESERVES_WHEN_RAMP_UP ? total_reserve[sc.name, gn, t] : 0.0) -
                prod_above[sc.name, gn, t-1] <= RU,
                base_name = "eq_ramp_up_$(sc.name)_$(gn)_$(t)"
        )

        eq_ramp_down[sc.name, gn, t] = @constraint(
                    model,
                    prod_above[sc.name, gn, t-1] +
                    (RESERVES_WHEN_RAMP_DOWN ? total_reserve[sc.name, gn, t-1] : 0.0) -
                    prod_above[sc.name, gn, t] <= RD,
                    base_name = "eq_ramp_dn_$(sc.name)_$(gn)_$(t)"
                )
    end
end


function _add_ramp_eqs_ori!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_prod_vars::Gar1962.ProdVars,
    formulation_ramping::MorLatRam2013.Ramping,
    formulation_status_vars::Gar1962.StatusVars,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_START_UP = true
    RESERVES_WHEN_RAMP_UP = true
    RESERVES_WHEN_RAMP_DOWN = true
    RESERVES_WHEN_SHUT_DOWN = true
    is_initially_on = (g.initial_status > 0)
    SU = g.startup_limit
    SD = g.shutdown_limit
    RU = g.ramp_up_limit
    RD = g.ramp_down_limit
    gn = g.name
    eq_ramp_down = _init(model, :eq_ramp_down)
    eq_ramp_up = _init(model, :eq_str_ramp_up)
    reserve = _total_reserves(model, g, sc, nCont, nInt)
    first_ramp_up = model[:first_ramp_up]
    first_ramp_dn = model[:first_ramp_dn]

    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    for t in 1:nCont+nInt
        time_invariant =
            (t > 1) ? (abs(g.min_power[t] - g.min_power[t-1]) < 1e-7) : true

        # Ramp up limit
        if t == 1
            if is_initially_on
                eq_ramp_up[sc.name, gn, t] = @constraint(
                    model,
                    g.min_power[t] +
                    prod_above[sc.name, gn, t] +
                    (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) <=
                    g.initial_power + RU
                    # , base_name="ori_rampup_$(gn)_$(t)"
                )
                first_ramp_up[gn] = eq_ramp_up[sc.name, gn, t]
            else
                first_ramp_up[gn] = @constraint(model, (g.min_power[t] + prod_above[sc.name, gn, t]) <= g.max_power[t])
            end
        else
            # amk: without accounting for time-varying min power terms,
            #      we might get an infeasible schedule, e.g. if min_power[t-1] = 0, min_power[t] = 10
            #      and ramp_up_limit = 5, the constraint (p'(t) + r(t) <= p'(t-1) + RU)
            #      would be satisfied with p'(t) = r(t) = p'(t-1) = 0
            #      Note that if switch_on[t] = 1, then eqns (20) or (21) go into effect
            if !time_invariant
                # Use equation (24) instead
                SU = g.startup_limit
                max_prod_this_period =
                    g.min_power[t] * is_on[gn, t] +
                    prod_above[sc.name, gn, t] +
                    (
                        RESERVES_WHEN_START_UP || RESERVES_WHEN_RAMP_UP ?
                        reserve[t] : 0.0
                    )
                min_prod_last_period =
                    g.min_power[t-1] * is_on[gn, t-1] +
                    prod_above[sc.name, gn, t-1]
                eq_ramp_up[gn, t] = @constraint(
                    model,
                    max_prod_this_period - min_prod_last_period <=
                    RU * is_on[gn, t-1] + SU * switch_on[gn, t]
                    # , base_name="ori_rampup_$(gn)_$(t)"
                )
            else
                # Equation (26) in Kneuven et al. (2020)
                # TODO: what if RU < SU? places too stringent upper bound
                # prod_above[gn, t] when starting up, and creates diff with (24).
                eq_ramp_up[sc.name, gn, t] = @constraint(
                    model,
                    prod_above[sc.name, gn, t] +
                    (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) -
                    prod_above[sc.name, gn, t-1] <= RU
                    # , base_name="ori_rampup_$(gn)_$(t)"
                )
            end
        end

        # Ramp down limit
        if t == 1
            if is_initially_on
                # TODO If RD < SD, or more specifically if
                #        min_power + RD < initial_power < SD
                #      then the generator should be able to shut down at time t = 1,
                #      but the constraint below will force the unit to produce power
                eq_ramp_down[sc.name, gn, t] = @constraint(
                    model,
                    g.initial_power -
                    # init_power_var[g.name] -
                    (g.min_power[t] + prod_above[sc.name, gn, t]) <= RD
                    # , base_name="ori_rampdn_$(gn)_$(t)"
                )
                first_ramp_dn[gn] = eq_ramp_down[sc.name, gn, t]
            else
                first_ramp_dn[gn] = @constraint(model, (g.min_power[t] + prod_above[sc.name, gn, t]) <= g.max_power[t])
            end

        else
            # amk: similar to ramp_up, need to account for time-dependent min_power
            if !time_invariant
                # Revert to (25)
                SD = g.shutdown_limit
                max_prod_last_period =
                    g.min_power[t-1] * is_on[gn, t-1] +
                    prod_above[sc.name, gn, t-1] +
                    (
                        RESERVES_WHEN_SHUT_DOWN || RESERVES_WHEN_RAMP_DOWN ?
                        reserve[t-1] : 0.0
                    )
                min_prod_this_period =
                    g.min_power[t] * is_on[gn, t] + prod_above[sc.name, gn, t]
                eq_ramp_down[gn, t] = @constraint(
                    model,
                    max_prod_last_period - min_prod_this_period <=
                    RD * is_on[gn, t] + SD * switch_off[gn, t]
                    # , base_name="ori_rampdn_$(gn)_$(t)"
                )
            else
                # Equation (27) in Kneuven et al. (2020)
                # TODO: Similar to above, what to do if shutting down in time t
                # and RD < SD? There is a difference with (25).
                eq_ramp_down[sc.name, gn, t] = @constraint(
                    model,
                    prod_above[sc.name, gn, t-1] +
                    (RESERVES_WHEN_RAMP_DOWN ? reserve[t-1] : 0.0) -
                    prod_above[sc.name, gn, t] <= RD
                    # , base_name="ori_rampdn_$(gn)_$(t)"
                )
            end
        end
    end
end



function _update_ramp_t0!(
    model::JuMP.Model,
    sc::ucRHScenario,
    init_power
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_START_UP = true
    RESERVES_WHEN_RAMP_UP = true
    RESERVES_WHEN_RAMP_DOWN = true
    RESERVES_WHEN_SHUT_DOWN = true
    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    # init_power_var = model[:init_power_var]

    fu = model[:first_ramp_up]
    fd = model[:first_ramp_dn]
    for g in sc.thermal_units
        gn = g.name
        RU = g.ramp_up_limit
        RD = g.ramp_down_limit
        SU = g.startup_limit
        SD = g.shutdown_limit
        reserve = _total_reserves_offset(model, g, sc, 1)

        t=1
        # Ramp up limit
        tmd2 = fu[gn]
        JuMP.delete(model,tmd2)
        fu[gn] = @constraint(
            model,
            g.min_power[t] +
            prod_above[sc.name, gn, t] +
            (RESERVES_WHEN_RAMP_UP ? reserve : 0.0) <=
            g.min_power[t] + init_power[gn] + RU
            # , base_name="rampup_$(gn)_$(t)"
        )

        # Ramp down limit
        tmd = fd[gn]
        JuMP.delete(model,tmd)
        fd[gn] = @constraint(
            model,
            (g.min_power[t] + init_power[gn]) -
            (g.min_power[t] + prod_above[sc.name, gn, t]) <= RD
            # , base_name="rampdn_$(gn)_$(t)"
        )
        # tmp = @constraint(model, prod_above[sc.name, gn, t] <= RD)
        # fd[gn] = tmp
    end
end

function _update_ramp_t0_var!(
    model::JuMP.Model,
    sc::ucRHScenario,
    init_power, offset
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_START_UP = true
    RESERVES_WHEN_RAMP_UP = true
    RESERVES_WHEN_RAMP_DOWN = true
    RESERVES_WHEN_SHUT_DOWN = true
    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    init_power_var = model[:init_power_var]
    init_power_var_dn = model[:init_power_var_dn]
    # reserve_var = model[:reserve_var]
    # reserve_var_dn = model[:reserve_var_dn]
    t = offset+1

    for g in sc.thermal_units
        gn = g.name
        RU = g.ramp_up_limit
        RD = g.ramp_down_limit
        SU = g.startup_limit
        SD = g.shutdown_limit
        reserve = _total_reserves_offset(model, g, sc, 1)
        # eq_ramp_up[sc.name, gn, t] = @constraint(
        #     model,
        #     g.min_power[t] +
        #     prod_above[sc.name, gn, t] +
        #     reserve_var[g.name,t] <=
        #     init_power_var[g.name] + RU
        # )

        JuMP.fix(init_power_var[gn], g.min_power[t]+init_power[gn],force=true)

        JuMP.fix(init_power_var_dn[gn], g.min_power[t]+init_power[gn],force=true)

        # tmp = @constraint(model, prod_above[sc.name, gn, t] <= RD)
        # fd[gn] = tmp
    end
end


function _update_ramp_t1!(
    model::JuMP.Model,
    g::ThermalUnit,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int,
    offset::Int
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_START_UP = true
    RESERVES_WHEN_RAMP_UP = true
    RESERVES_WHEN_RAMP_DOWN = true
    RESERVES_WHEN_SHUT_DOWN = true
    is_initially_on = (g.initial_status > 0)
    SU = g.startup_limit
    SD = g.shutdown_limit
    RU = g.ramp_up_limit
    RD = g.ramp_down_limit
    gn = g.name
    # eq_ramp_down = _init(model, :eq_ramp_down)
    # eq_ramp_up = _init(model, :eq_str_ramp_up)
    reserve = _total_reserves_offset2(model, g, sc, nCont, nInt, offset)
    eq_ramp_down = model[:eq_ramp_down]
    eq_ramp_up = model[:eq_str_ramp_up]

    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    for t in 2:nCont+nInt
        time_invariant =
            (t > 1) ? (abs(g.min_power[t] - g.min_power[t-1]) < 1e-7) : true

        # Ramp up limit
        # remove constraint
        # amk: without accounting for time-varying min power terms,
        #      we might get an infeasible schedule, e.g. if min_power[t-1] = 0, min_power[t] = 10
        #      and ramp_up_limit = 5, the constraint (p'(t) + r(t) <= p'(t-1) + RU)
        #      would be satisfied with p'(t) = r(t) = p'(t-1) = 0
        #      Note that if switch_on[t] = 1, then eqns (20) or (21) go into effect
        if !time_invariant
            # Use equation (24) instead
            SU = g.startup_limit
            max_prod_this_period =
                g.min_power[t] * is_on[gn, t] +
                prod_above[sc.name, gn, t] +
                (
                    RESERVES_WHEN_START_UP || RESERVES_WHEN_RAMP_UP ?
                    reserve[t] : 0.0
                )
            min_prod_last_period =
                g.min_power[t-1] * is_on[gn, t-1] +
                prod_above[sc.name, gn, t-1]
            tmd = eq_ramp_up[gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_up[gn, t] = @constraint(
                model,
                max_prod_this_period - min_prod_last_period <=
                RU * is_on[gn, t-1] + SU * switch_on[gn, t]
                # , base_name="rampup_$(gn)_$(t)"
            )
        else
            # Equation (26) in Kneuven et al. (2020)
            # TODO: what if RU < SU? places too stringent upper bound
            # prod_above[gn, t] when starting up, and creates diff with (24).
            tmd = eq_ramp_up[sc.name, gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_up[sc.name, gn, t] = @constraint(
                model,
                prod_above[sc.name, gn, t] +
                (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) -
                prod_above[sc.name, gn, t-1] <= RU
                # , base_name="rampup_$(gn)_$(t)"
            )
        end

        # Ramp down limit
        # amk: similar to ramp_up, need to account for time-dependent min_power
        if !time_invariant
            # Revert to (25)
            SD = g.shutdown_limit
            max_prod_last_period =
                g.min_power[t-1] * is_on[gn, t-1] +
                prod_above[sc.name, gn, t-1] +
                (
                    RESERVES_WHEN_SHUT_DOWN || RESERVES_WHEN_RAMP_DOWN ?
                    reserve[t-1] : 0.0
                )
            min_prod_this_period =
                g.min_power[t] * is_on[gn, t] + prod_above[sc.name, gn, t]
            tmd = eq_ramp_down[gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_down[gn, t] = @constraint(
                model,
                max_prod_last_period - min_prod_this_period <=
                RD * is_on[gn, t] + SD * switch_off[gn, t]
                # , base_name="rampdn_$(gn)_$(t)"
            )
        else
            # Equation (27) in Kneuven et al. (2020)
            # TODO: Similar to above, what to do if shutting down in time t
            # and RD < SD? There is a difference with (25).
            tmd = eq_ramp_down[sc.name, gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_down[sc.name, gn, t] = @constraint(
                model,
                prod_above[sc.name, gn, t-1] +
                (RESERVES_WHEN_RAMP_DOWN ? reserve[t-1] : 0.0) -
                prod_above[sc.name, gn, t] <= RD
                # , base_name="rampdn_$(gn)_$(t)"
            )
        end
    end
end


function _update_ramp_t1_var!(
    model::JuMP.Model,
    g::ThermalUnit,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int,
    offset::Int
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_START_UP = true
    RESERVES_WHEN_RAMP_UP = true
    RESERVES_WHEN_RAMP_DOWN = true
    RESERVES_WHEN_SHUT_DOWN = true
    is_initially_on = (g.initial_status > 0)
    SU = g.startup_limit
    SD = g.shutdown_limit
    RU = g.ramp_up_limit
    RD = g.ramp_down_limit
    gn = g.name
    # eq_ramp_down = _init(model, :eq_ramp_down)
    # eq_ramp_up = _init(model, :eq_str_ramp_up)
    reserve = _total_reserves_offset2(model, g, sc, nCont, nInt, offset)
    eq_ramp_down = model[:eq_ramp_down]
    eq_ramp_up = model[:eq_str_ramp_up]

    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    for t in 2:nCont+nInt
        time_invariant =
            (t > 1) ? (abs(g.min_power[t] - g.min_power[t-1]) < 1e-7) : true

        # Ramp up limit
        # remove constraint
        # amk: without accounting for time-varying min power terms,
        #      we might get an infeasible schedule, e.g. if min_power[t-1] = 0, min_power[t] = 10
        #      and ramp_up_limit = 5, the constraint (p'(t) + r(t) <= p'(t-1) + RU)
        #      would be satisfied with p'(t) = r(t) = p'(t-1) = 0
        #      Note that if switch_on[t] = 1, then eqns (20) or (21) go into effect
        if !time_invariant
            # Use equation (24) instead
            SU = g.startup_limit
            max_prod_this_period =
                g.min_power[t] * is_on[gn, t] +
                prod_above[sc.name, gn, t] +
                (
                    RESERVES_WHEN_START_UP || RESERVES_WHEN_RAMP_UP ?
                    reserve[t] : 0.0
                )
            min_prod_last_period =
                g.min_power[t-1] * is_on[gn, t-1] +
                prod_above[sc.name, gn, t-1]
            tmd = eq_ramp_up[gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_up[gn, t] = @constraint(
                model,
                max_prod_this_period - min_prod_last_period <=
                RU * is_on[gn, t-1] + SU * switch_on[gn, t]
                # , base_name="rampup_$(gn)_$(t)"
            )
        else
            # Equation (26) in Kneuven et al. (2020)
            # TODO: what if RU < SU? places too stringent upper bound
            # prod_above[gn, t] when starting up, and creates diff with (24).
            tmd = eq_ramp_up[sc.name, gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_up[sc.name, gn, t] = @constraint(
                model,
                prod_above[sc.name, gn, t] +
                (RESERVES_WHEN_RAMP_UP ? reserve[t] : 0.0) -
                prod_above[sc.name, gn, t-1] <= RU
                # , base_name="rampup_$(gn)_$(t)"
            )
        end

        # Ramp down limit
        # amk: similar to ramp_up, need to account for time-dependent min_power
        if !time_invariant
            # Revert to (25)
            SD = g.shutdown_limit
            max_prod_last_period =
                g.min_power[t-1] * is_on[gn, t-1] +
                prod_above[sc.name, gn, t-1] +
                (
                    RESERVES_WHEN_SHUT_DOWN || RESERVES_WHEN_RAMP_DOWN ?
                    reserve[t-1] : 0.0
                )
            min_prod_this_period =
                g.min_power[t] * is_on[gn, t] + prod_above[sc.name, gn, t]
            tmd = eq_ramp_down[gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_down[gn, t] = @constraint(
                model,
                max_prod_last_period - min_prod_this_period <=
                RD * is_on[gn, t] + SD * switch_off[gn, t]
                # , base_name="rampdn_$(gn)_$(t)"
            )
        else
            # Equation (27) in Kneuven et al. (2020)
            # TODO: Similar to above, what to do if shutting down in time t
            # and RD < SD? There is a difference with (25).
            tmd = eq_ramp_down[sc.name, gn, t]
            JuMP.delete(model,tmd)
            eq_ramp_down[sc.name, gn, t] = @constraint(
                model,
                prod_above[sc.name, gn, t-1] +
                (RESERVES_WHEN_RAMP_DOWN ? reserve[t-1] : 0.0) -
                prod_above[sc.name, gn, t] <= RD
                # , base_name="rampdn_$(gn)_$(t)"
            )
        end
    end
end