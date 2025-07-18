# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

function _add_ramp_eqs!(
    model::JuMP.Model,
    g::ThermalUnit,
    formulation_prod_vars::Gar1962.ProdVars,
    formulation_ramping::PanGua2016.Ramping,
    formulation_status_vars::Gar1962.StatusVars,
    sc::ucRHScenario,
    nCont::Int,
    nInt::Int
)::Nothing
    # TODO: Move upper case constants to model[:instance]
    RESERVES_WHEN_SHUT_DOWN = true
    gn = g.name
    reserve = _total_reserves(model, g, sc, nCont, nInt)
    eq_str_prod_limit = _init(model, :eq_str_prod_limit)
    eq_prod_limit_ramp_up_extra_period =
        _init(model, :eq_prod_limit_ramp_up_extra_period)
    eq_prod_limit_shutdown_trajectory =
        _init(model, :eq_prod_limit_shutdown_trajectory)
    UT = g.min_uptime
    SU = g.startup_limit   # startup rate, i.e., max production right after startup
    SD = g.shutdown_limit  # shutdown rate, i.e., max production right before shutdown
    RU = g.ramp_up_limit   # ramp up rate
    RD = g.ramp_down_limit # ramp down rate
    T = nCont+nInt

    # Gar1962.ProdVars
    prod_above = model[:prod_above]

    # Gar1962.StatusVars
    is_on = model[:is_on]
    switch_off = model[:switch_off]
    switch_on = model[:switch_on]

    for t in 1:T
        Pbar = g.max_power[t]
        if Pbar < 1e-7
            # Skip this time period if max power = 0
            continue
        end

        #TRD = floor((Pbar - SU) / RD) # ramp down time
        # TODO check amk changed TRD wrt Kneuven et al.
        TRD = ceil((Pbar - SD) / RD)  # ramp down time
        TRU = floor((Pbar - SU) / RU) # ramp up time, can be negative if Pbar < SU

        # TODO check initial time periods: what if generator has been running for x periods?
        # But maybe ok as long as (35) and (36) are also used...
        if UT > 1
            # Equation (38) in Kneuven et al. (2020)
            # Generalization of (20)
            # Necessary that if any of the switch_on = 1 in the sum,
            # then switch_off[gn, t+1] = 0
            eq_str_prod_limit[sc.name, gn, t] = @constraint(
                model,
                prod_above[sc.name, gn, t] +
                g.min_power[t] * is_on[gn, t] +
                reserve[t] <=
                Pbar * is_on[gn, t] -
                (t < T ? (Pbar - SD) * switch_off[gn, t+1] : 0.0) - sum(
                    (Pbar - (SU + i * RU)) * switch_on[gn, t-i] for
                    i in 0:min(UT - 2, TRU, t - 1)
                )
                # ,base_name = "ramp_Kneuven_$(gn)_$(t)"
            )

            if UT - 2 < TRU
                # Equation (40) in Kneuven et al. (2020)
                # Covers an additional time period of the ramp-up trajectory, compared to (38)
                eq_prod_limit_ramp_up_extra_period[sc.name, gn, t] =
                    @constraint(
                        model,
                        prod_above[sc.name, gn, t] +
                        g.min_power[t] * is_on[gn, t] +
                        reserve[t] <=
                        Pbar * is_on[gn, t] - sum(
                            (Pbar - (SU + i * RU)) * switch_on[gn, t-i] for
                            i in 0:min(UT - 1, TRU, t - 1)
                        )
                        # ,base_name = "ramp_Kneuven_trup_$(gn)_$(t)"
                    )
            end

            # Add in shutdown trajectory if KSD >= 0 (else this is dominated by (38))
            KSD = min(TRD, UT - 1, T - t - 1)
            if KSD > 0
                KSU = min(TRU, UT - 2 - KSD, t - 1)
                # Equation (41) in Kneuven et al. (2020)
                eq_prod_limit_shutdown_trajectory[sc.name, gn, t] = @constraint(
                    model,
                    prod_above[sc.name, gn, t] +
                    g.min_power[t] * is_on[gn, t] +
                    (RESERVES_WHEN_SHUT_DOWN ? reserve[t] : 0.0) <=
                    Pbar * is_on[gn, t] - sum(
                        (Pbar - (SD + i * RD)) * switch_off[gn, t+1+i] for
                        i in 0:KSD
                    ) - sum(
                        (Pbar - (SU + i * RU)) * switch_on[gn, t-i] for
                        i in 0:KSU
                    ) - (
                        (KSU >= TRU || KSU > t - 2) ? 0.0 :
                        max(0, (SU + (KSU + 1) * RU) - (SD + TRD * RD)) *
                        switch_on[gn, t-(KSU+1)]
                    )
                    # ,base_name = "ramp_Kneuven_trdn_$(gn)_$(t)"
                )
            end
        end
    end
end
