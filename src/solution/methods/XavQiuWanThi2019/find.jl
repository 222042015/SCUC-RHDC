# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

 


import Base.Threads: @threads

function _find_violations(
    model::JuMP.Model,
    sc::ucRHScenario;
    max_per_line::Int,
    max_per_period::Int,
    nStart::Int,
    nEnd::Int
)
    @info "called this one"
    instance = model[:instance]
    net_injection = model[:net_injection]
    overflow = model[:overflow]
    length(sc.buses) > 1 || return []
    violations = []

    non_slack_buses = [b for b in sc.buses if b.offset > 0]
    net_injection_values = [
        value(net_injection[sc.name, b.name, t]) for b in non_slack_buses,
        t in nStart:nEnd
    ]
    overflow_values = [
        value(overflow[sc.name, lm.name, t]) for lm in sc.lines,
        t in nStart:nEnd
    ]
    violations = ucRH._find_violations(
        instance = instance,
        sc = sc,
        net_injections = net_injection_values,
        overflow = overflow_values,
        isf = sc.isf,
        lodf = sc.lodf,
        max_per_line = max_per_line,
        max_per_period = max_per_period,
        nStart=nStart,
        nEnd=nEnd,
        offset=model[:offset],
    )
    return violations
end


function _find_violations_callback(
    model::JuMP.Model,
    sc::ucRHScenario;
    max_per_line::Int,
    max_per_period::Int,
    nStart::Int,
    nEnd::Int,
    cb_data
)
    @info "called this one"
    instance = model[:instance]
    net_injection = model[:net_injection]
    overflow = model[:overflow]
    length(sc.buses) > 1 || return []
    violations = []
    non_slack_buses = [b for b in sc.buses if b.offset > 0]
    net_injection_values = [
        callback_value(cb_data,net_injection[sc.name, b.name, t]) for b in non_slack_buses,
        t in nStart:nEnd
    ]
    overflow_values = [
        callback_value(cb_data,overflow[sc.name, lm.name, t]) for lm in sc.lines,
        t in nStart:nEnd
    ]
    violations = ucRH._find_violations(
        instance = instance,
        sc = sc,
        net_injections = net_injection_values,
        overflow = overflow_values,
        isf = sc.isf,
        lodf = sc.lodf,
        max_per_line = max_per_line,
        max_per_period = max_per_period,
        nStart=nStart,
        nEnd=nEnd,
        offset=model[:offset]
    )
    return violations
end

"""
    function _find_violations(
        instance::ucRHInstance,
        net_injections::Array{Float64, 2};
        isf::Array{Float64,2},
        lodf::Array{Float64,2},
        max_per_line::Int,
        max_per_period::Int,
    )::Array{_Violation, 1}

Find transmission constraint violations (both pre-contingency, as well as
post-contingency).

The argument `net_injection` should be a (B-1) x T matrix, where B is the
number of buses and T is the number of time periods. The arguments `isf` and
`lodf` can be computed using ucRH.injection_shift_factors and
ucRH.line_outage_factors. The argument `overflow` specifies how much
flow above the transmission limits (in MW) is allowed. It should be an L x T
matrix, where L is the number of transmission lines.
"""
function _find_violations(;
    instance::ucRHInstance,
    sc::ucRHScenario,
    net_injections::Array{Float64,2},
    overflow::Array{Float64,2},
    isf::Array{Float64,2},
    lodf::Array{Float64,2},
    max_per_line::Int,
    max_per_period::Int,
    nStart::Int,
    nEnd::Int,
    offset::Int
)::Array{_Violation,1}
    B = length(sc.buses) - 1
    L = length(sc.lines)
    T = nEnd - nStart+1
    printstyled("Calling Find Violation with T=$(T), ($(nStart)~$(nEnd))\n"; color=:red)
    @info "Calling Find Violation with T=$(T), ($(nStart)~$(nEnd))\n"
    K = nthreads()
    println("number of threads: $(K)")

    size(net_injections) == (B, T) || error("net_injections has incorrect size")
    size(isf) == (L, B) || error("isf has incorrect size")
    size(lodf) == (L, L) || error("lodf has incorrect size")

    filters = Dict(
        t => _ViolationFilter(
            max_total = max_per_period,
            max_per_line = max_per_line,
        ) for t in 1:T
    )

    pre_flow::Array{Float64} = zeros(L, K)           # pre_flow[lm, thread]
    # post_flow::Array{Float64} = zeros(L, L, K)       # post_flow[lm, lc, thread]
    pre_v::Array{Float64} = zeros(L, K)              # pre_v[lm, thread]
    # post_v::Array{Float64} = zeros(L, L, K)          # post_v[lm, lc, thread]

    normal_limits::Array{Float64,2} = [
        l.normal_flow_limit[t+offset] + overflow[l.offset, t] for l in sc.lines,
        t in 1:T
    ]

    emergency_limits::Array{Float64,2} = [
        l.emergency_flow_limit[t+offset] + overflow[l.offset, t] for l in sc.lines,
        t in 1:T
    ]

    is_vulnerable::Array{Bool} = zeros(Bool, L)
    for c in sc.contingencies
        is_vulnerable[c.lines[1].offset] = true
    end

    @threads for t in 1:T
        k = threadid()

        # Pre-contingency flows
        pre_flow[:, k] = isf * net_injections[:, t]

        # Post-contingency flows
        # for lc in 1:L, lm in 1:L
        #     post_flow[lm, lc, k] =
        #         pre_flow[lm, k] + pre_flow[lc, k] * lodf[lm, lc]
        # end

        to_post = []

        t1 = @elapsed begin
            # Pre-contingency violations
            for lm in 1:L
                # pre_v[lm, k] = max(
                #     0.0,
                #     pre_flow[lm, k] - normal_limits[lm, t],
                #     -pre_flow[lm, k] - normal_limits[lm, t],
                # )

                tmp = pre_flow[lm, k]
                if tmp < 0.0
                    tmp = -tmp
                end
                tmp -= normal_limits[lm, t]
                if tmp < 0.0
                    tmp = 0.0
                end
                # TODO:: skip 1e-5

                pre_v[lm, k] = tmp
            end

            # Post-contingency violations
            for lc in 1:L
                if !is_vulnerable[lc]
                    continue
                end

                tmp2 = pre_flow[lc, k]
                for lm in 1:L
                    # post_v[lm, lc, k] = max(
                    #     0.0,
                    #     post_flow[lm, lc, k] - emergency_limits[lm, t],
                    #     -post_flow[lm, lc, k] - emergency_limits[lm, t],
                    # )

                    tmp = pre_flow[lm, k] + tmp2 * lodf[lm, lc]
                    if tmp < 0.0
                        tmp = -tmp
                    end
                    tmp -= emergency_limits[lm, t]

                    if tmp > 1e-5 
                        push!(to_post, [[lm, lc], tmp])
                    end
                end
            end
        end

        t2 = @elapsed begin
            # Offer pre-contingency violations
            for lm in 1:L
                if pre_v[lm, k] > 1e-5
                    _offer(
                        filters[t],
                        _Violation(
                            time = t,
                            monitored_line = sc.lines[lm],
                            outage_line = nothing,
                            amount = pre_v[lm, k],
                        ),
                    )
                end
            end
        end

        t3 = @elapsed begin
            # Offer post-contingency violations
            # for lm in 1:L, lc in 1:L
            #     if post_v[lm, lc, k] > 1e-5 && is_vulnerable[lc]
            #         _offer(
            #             filters[t],
            #             _Violation(
            #                 time = t,
            #                 monitored_line = sc.lines[lm],
            #                 outage_line = sc.lines[lc],
            #                 amount = post_v[lm, lc, k],
            #             ),
            #         )
            #     end
            # end
            for ele in to_post
                lm = ele[1][1]
                lc = ele[1][2]
                vv = ele[2]
                # println("$(lm) $(lc) $(vv)")
                _offer(
                        filters[t],
                        _Violation(
                            time = t,
                            monitored_line = sc.lines[lm],
                            outage_line = sc.lines[lc],
                            amount = vv,
                            # amount = post_v[lm, lc, k],
                        ),
                    )
            end
        end
        # println("!!!!!!!!!!!!! 3 times:::::: $(t1) $(t2) $(t3)")
    end

    violations = _Violation[]
    for t in 1:T
        append!(violations, _query(filters[t]))
    end

    return violations
end



function _find_violations_ori(;
    instance::ucRHInstance,
    sc::ucRHScenario,
    net_injections::Array{Float64,2},
    overflow::Array{Float64,2},
    isf::Array{Float64,2},
    lodf::Array{Float64,2},
    max_per_line::Int,
    max_per_period::Int,
)::Array{_Violation,1}
    B = length(sc.buses) - 1
    L = length(sc.lines)
    T = instance.time
    K = nthreads()

    size(net_injections) == (B, T) || error("net_injections has incorrect size")
    size(isf) == (L, B) || error("isf has incorrect size")
    size(lodf) == (L, L) || error("lodf has incorrect size")

    filters = Dict(
        t => _ViolationFilter(
            max_total = max_per_period,
            max_per_line = max_per_line,
        ) for t in 1:T
    )

    pre_flow::Array{Float64} = zeros(L, K)           # pre_flow[lm, thread]
    post_flow::Array{Float64} = zeros(L, L, K)       # post_flow[lm, lc, thread]
    pre_v::Array{Float64} = zeros(L, K)              # pre_v[lm, thread]
    post_v::Array{Float64} = zeros(L, L, K)          # post_v[lm, lc, thread]

    normal_limits::Array{Float64,2} = [
        l.normal_flow_limit[t] + overflow[l.offset, t] for l in sc.lines,
        t in 1:T
    ]

    emergency_limits::Array{Float64,2} = [
        l.emergency_flow_limit[t] + overflow[l.offset, t] for l in sc.lines,
        t in 1:T
    ]

    is_vulnerable::Array{Bool} = zeros(Bool, L)
    for c in sc.contingencies
        is_vulnerable[c.lines[1].offset] = true
    end

    @threads for t in 1:T
        k = threadid()

        # Pre-contingency flows
        pre_flow[:, k] = isf * net_injections[:, t]

        # Post-contingency flows
        for lc in 1:L, lm in 1:L
            post_flow[lm, lc, k] =
                pre_flow[lm, k] + pre_flow[lc, k] * lodf[lm, lc]
        end

        # Pre-contingency violations
        for lm in 1:L
            pre_v[lm, k] = max(
                0.0,
                pre_flow[lm, k] - normal_limits[lm, t],
                -pre_flow[lm, k] - normal_limits[lm, t],
            )
        end

        # Post-contingency violations
        for lc in 1:L, lm in 1:L
            post_v[lm, lc, k] = max(
                0.0,
                post_flow[lm, lc, k] - emergency_limits[lm, t],
                -post_flow[lm, lc, k] - emergency_limits[lm, t],
            )
        end

        # Offer pre-contingency violations
        for lm in 1:L
            if pre_v[lm, k] > 1e-5
                _offer(
                    filters[t],
                    _Violation(
                        time = t,
                        monitored_line = sc.lines[lm],
                        outage_line = nothing,
                        amount = pre_v[lm, k],
                    ),
                )
            end
        end

        # Offer post-contingency violations
        for lm in 1:L, lc in 1:L
            if post_v[lm, lc, k] > 1e-5 && is_vulnerable[lc]
                _offer(
                    filters[t],
                    _Violation(
                        time = t,
                        monitored_line = sc.lines[lm],
                        outage_line = sc.lines[lc],
                        amount = post_v[lm, lc, k],
                    ),
                )
            end
        end
    end

    violations = _Violation[]
    for t in 1:instance.time
        append!(violations, _query(filters[t]))
    end

    return violations
end
