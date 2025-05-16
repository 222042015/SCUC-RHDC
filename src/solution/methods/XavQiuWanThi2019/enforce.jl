# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

function _enforce_transmission(
    model::JuMP.Model,
    violations::Vector{_Violation},
    sc::ucRHScenario,
)::Nothing
    for v in violations
        _enforce_transmission(
            model = model,
            sc = sc,
            violation = v,
            isf = sc.isf,
            lodf = sc.lodf,
        )
    end
    return
end

function _enforce_transmission(;
    model::JuMP.Model,
    sc::ucRHScenario,
    violation::_Violation,
    isf::Matrix{Float64},
    lodf::Matrix{Float64},
)::Nothing
    instance = model[:instance]
    limit::Float64 = 0.0
    overflow = model[:overflow]
    net_injection = model[:net_injection]

    if violation.outage_line === nothing
        limit = violation.monitored_line.normal_flow_limit[violation.time]
        @info @sprintf(
            "    %8.3f MW overflow in %-5s time %3d (pre-contingency, scenario %s)",
            violation.amount,
            violation.monitored_line.name,
            violation.time,
            sc.name,
        )
    else
        limit = violation.monitored_line.emergency_flow_limit[violation.time]
        @info @sprintf(
            "    %8.3f MW overflow in %-5s time %3d (outage: line %s, scenario %s)",
            violation.amount,
            violation.monitored_line.name,
            violation.time,
            violation.outage_line.name,
            sc.name,
        )
    end

    fm = violation.monitored_line.name
    t = violation.time
    flow = @variable(model, base_name = "flow[$fm,$t]")

    v = overflow[sc.name, violation.monitored_line.name, violation.time]
    @constraint(model, flow <= limit + v)
    @constraint(model, -flow <= limit + v)

    if violation.outage_line === nothing
        @constraint(
            model,
            flow == sum(
                net_injection[sc.name, b.name, violation.time] *
                isf[violation.monitored_line.offset, b.offset] for
                b in sc.buses if b.offset > 0
            )
        )
    else
        @constraint(
            model,
            flow == sum(
                net_injection[sc.name, b.name, violation.time] * (
                    isf[violation.monitored_line.offset, b.offset] + (
                        lodf[
                            violation.monitored_line.offset,
                            violation.outage_line.offset,
                        ] * isf[violation.outage_line.offset, b.offset]
                    )
                ) for b in sc.buses if b.offset > 0
            )
        )
    end
    return nothing
end


function _enforce_transmission_callback(
    model::JuMP.Model,
    violations::Vector{_Violation},
    sc::ucRHScenario,
    cb_data
)::Nothing
    for v in violations
        _enforce_transmission_callback(
            model = model,
            sc = sc,
            violation = v,
            isf = sc.isf,
            lodf = sc.lodf,
            cb_data = cb_data
        )
    end
    return
end

function _enforce_transmission_callback(;
    model::JuMP.Model,
    sc::ucRHScenario,
    violation::_Violation,
    isf::Matrix{Float64},
    lodf::Matrix{Float64},
    cb_data
)::Nothing
    use_restore = false

    instance = model[:instance]
    limit::Float64 = 0.0
    overflow = model[:overflow]
    net_injection = model[:net_injection]

    if violation.outage_line === nothing
        limit = violation.monitored_line.normal_flow_limit[violation.time]
        @info @sprintf(
            "    %8.3f MW CBBB overflow in %-5s time %3d (pre-contingency, scenario %s)",
            violation.amount,
            violation.monitored_line.name,
            violation.time,
            sc.name,
        )
    else
        limit = violation.monitored_line.emergency_flow_limit[violation.time]
        @info @sprintf(
            "    %8.3f MW CBBB  overflow in %-5s time %3d (outage: line %s, scenario %s)",
            violation.amount,
            violation.monitored_line.name,
            violation.time,
            violation.outage_line.name,
            sc.name,
        )
    end

    fm = violation.monitored_line.name
    t = violation.time
    # flow = @variable(model, base_name = "flow[$fm,$t]")

    v = overflow[sc.name, violation.monitored_line.name, violation.time]
    # tc = @constraint(model, flow <= limit + v)
    # tc = @constraint(model, -flow <= limit + v)

    if violation.outage_line === nothing
        con = @build_constraint(sum(
                net_injection[sc.name, b.name, violation.time] *
                isf[violation.monitored_line.offset, b.offset] for
                b in sc.buses if b.offset > 0
            ) <= limit + v
        )
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)

        con = @build_constraint(sum(
                net_injection[sc.name, b.name, violation.time] *
                isf[violation.monitored_line.offset, b.offset] for
                b in sc.buses if b.offset > 0
            ) >= -limit - v
        )
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)

        # tc = @constraint(
        #     model,
        #     flow == sum(
        #         net_injection[sc.name, b.name, violation.time] *
        #         isf[violation.monitored_line.offset, b.offset] for
        #         b in sc.buses if b.offset > 0
        #     )
        # )
    else
        tmp_exp = AffExpr()
        tmp_ex = []
        for b in sc.buses 
            if b.offset > 0
                a = lodf[violation.monitored_line.offset,violation.outage_line.offset,]
                bz = isf[violation.monitored_line.offset, b.offset]
                c = isf[violation.outage_line.offset, b.offset]
                ed = bz+a*c
                if (abs(ed)>1e-6)
                    add_to_expression!(tmp_exp, ed, net_injection[sc.name, b.name, violation.time])
                    if use_restore
                        push!(tmp_ex, [ed,"netj",sc.name, b.name, violation.time])
                    end
                end
            end
        end

        con = @build_constraint(tmp_exp >= -limit - v)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        con = @build_constraint(tmp_exp <= limit + v)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)

        # tc = @constraint(
        #     model,
        #     flow == tmp_exp
        # )
    end
    return nothing
end
