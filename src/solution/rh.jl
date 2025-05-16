
# function rh_full!(
#     model,
#     binStart,
#     binEnd,
#     contEnd,
#     to_fix
# )
#     # 

# end



function run_rh!(
    model,
    nCont,
    nInt,
    offset
)
    JuMP.set_optimizer_attribute(model, "Threads", 1)
    to_fix = []
    # get global info and variables
    instance = model[:instance].scenarios[1]
    variables = JuMP.all_variables(model)

    # optimize the bin-cont model
    # MOI.set(model, MOI.NumberOfThreads(), 1)
    # JuMP.set_time_limit_sec(model, 1000.0)
    adc = _init(model, :adc)
    adv = _init(model, :adv)
    ucRH.optimize!(model,1,nCont+nInt)

    ##########################################################################
    T = nCont+nInt
    up_reserved = []
    down_reserved = []
    last_stat = []
    p_last = Dict()
    # record solution
    for b in instance.thermal_units
        # record global solution 
        mup = b.min_uptime
        mdn = b.min_downtime
        # @info "$(b.name) min up time: $(mup)    min down time: $(mdn)"
        minup_flag=false
        mindn_flag=false

        
        push!(last_stat, ["is_on", b.name, nInt+offset, value(model[:is_on][b.name, nInt])])
        for t in 1:nInt
            vfix = value(model[:is_on][b.name, t])
            push!(to_fix, ["is_on", b.name, t+offset, vfix])
            vswon = value(model[:switch_on][b.name, t])
            # push!(to_fix, ["switch_on", b.name, t+offset, vswon])
            vswof = value(model[:switch_off][b.name, t])
            # push!(to_fix, ["switch_off", b.name, t+offset, vswof])

            # if "g56"==b.name
            #     printstyled("$(b.name)_$(t): $(vswon) $(vswof) $(vfix)\n"; color=:yellow)
            # end

            # check if need to update reserved down
            if !minup_flag && vswon>0.5 && t+mup-1>nInt
                # need to add to add to up_reserved
                to_add = t+mup-1-T
                push!(up_reserved,[b.name, to_add])
                minup_flag=true
                # @info "---At t=$(t), $(b.name) turned on, with minup=$(mup), need to be on for $(to_add) more days after T=$(nInt)"
            end
            if !mindn_flag && vswof>0.5 && t+mdn-1>nInt 
                # need to add to add to up_reserved
                to_add = t+mdn-1-T
                push!(down_reserved,[b.name, to_add])
                mindn_flag=true
                # @info "---At t=$(t), $(b.name) turned off, with minup=$(mdn), need to be down for $(to_add) more days after T=$(nInt)"
            end

        end
        # store the power output 
        vpabov = value(model[:prod_above][instance.name, b.name, nInt])
        # push!(p_last, [b.name, vpabov])
        p_last[b.name] = vpabov
        # @info "--- $(b.name)'s power output at t=$(nInt) is $(vpabov)"


        # for t in nInt+1:nCont+nInt
        #     vfix = value(model[:is_on][b.name, t])
        #     vswon = value(model[:switch_on][b.name, t])
        #     vswof = value(model[:switch_off][b.name, t])

        #     # check if need to update reserved down

        #     # check if need to update reserved down
        #     if !minup_flag && vswon>0.5 && t+mup-1>T 
        #         # need to add to add to up_reserved
        #         to_add = t+mup-1-T
        #         push!(up_reserved,[b.name, to_add])
        #         minup_flag=true
        #         @info "---At t=$(t), $(b.name) turned on, with minup=$(mup), need to be on for $(to_add) more days after T=$(T)"
        #     end
        #     if !mindn_flag && vswof>0.5 && t+mdn-1>T 
        #         # need to add to add to up_reserved
        #         to_add = t+mdn-1-T
        #         push!(down_reserved,[b.name, to_add])
        #         mindn_flag=true
        #         @info "---At t=$(t), $(b.name) turned off, with minup=$(mdn), need to be down for $(to_add) more days after T=$(T)"
        #     end
        # end
    end

    # TODO:: add update so the window moves forward

    return to_fix, up_reserved, down_reserved, p_last, last_stat
end


# function rh_full!(
#     model,
#     binStart,
#     binEnd,
#     contEnd,
#     to_fix
# )
#     # 

# end

# julia --threads 16 test_cb.jl | tee /home/lxyang/SCUC/julia/package/gittemp/ucRH/logs/case6468rte_cb2_16thread.log
# julia --threads 8 test_cb.jl | tee /home/lxyang/SCUC/julia/package/gittemp/ucRH/logs/case6468rte_cb2_1thread.log

function run_rh_cb!(
    model,
    nCont,
    nInt,
    stepsize,
    offset
)   
    fix_startup = false
    # get global info and variables
    instance = model[:instance].scenarios[1]
    is_on_ts = model[:is_on_ts]
    switch_on_ts = model[:switch_on_ts]
    switch_off_ts = model[:switch_off_ts]
    startup_ts = model[:startup_ts]
    initial_total_reserve = model[:initial_total_reserve]

    ##########################################################################
    # T = nCont+nInt
    if model[:rh]
        JuMP.optimize!(model)
    else
        ucRH.optimize!(model, nStart=1, nEnd=nCont+nInt)
    end
    # JuMP.write_to_file(model, "rh_cb.lp")

    # @info "wrote wrong lp"
    # record solution
    for g in instance.thermal_units # iterate over all thermal units
        init_status = g.initial_status
        for t in 1:stepsize # only the solutions from [offset, offset+stepsize] are set to fixed in the next step
            ison = value(model[:is_on][g.name, t])
            is_on_ts[g.name, offset+t] = ison
            switch_on_ts[g.name, offset+t] = value(model[:switch_on][g.name, t])
            switch_off_ts[g.name, offset+t] = value(model[:switch_off][g.name, t])
            
            # update the global stauts of g, initial_status.
            # integer, > 0, number of time steps g is on before current time step
            # integer, < 0, number of time steps g is off before current time step
            if init_status > 0
                if ison > 0.5
                    init_status = init_status + 1
                else
                    init_status = -1
                end
            else
                if ison > 0.5
                    init_status = 1
                else
                    init_status = init_status - 1
                end
                
            end

            if fix_startup
                S = length(g.startup_categories)
                for s in 1:S
                    startup_ts[g.name, offset+t, s] = value(model[:startup][g.name, t, s])
                end
            end
        end

        

        # update the instance inplace
        g.initial_status = init_status
        g.initial_power = value(model[:prod_above][instance.name, g.name, stepsize] + g.min_power[offset+stepsize] * value(model[:is_on][g.name, stepsize]))
        spinning_reserves = [r for r in g.reserves if r.type == "spinning"]
        if !isempty(spinning_reserves)
            initial_total_reserve[g.name] = sum(value(model[:reserve][instance.name, r.name, g.name, stepsize]) for r in spinning_reserves)
        else
            initial_total_reserve[g.name] = 0.0
        end
    end
end