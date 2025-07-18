# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020-2021, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

module XavQiuAhm2021

using Distributions
import ..ucRHInstance
import ..ucRHScenario

"""
    struct Randomization
        cost = Uniform(0.95, 1.05)
        load_profile_mu = [...]
        load_profile_sigma = [...]
        load_share = Uniform(0.90, 1.10)
        peak_load = Uniform(0.6 * 0.925, 0.6 * 1.075)
        randomize_costs = true
        randomize_load_profile = true
        randomize_load_share = true
    end

Randomization method that changes: (1) production and startup costs, (2)
share of load coming from each bus, (3) peak system load, and (4) temporal
load profile, as follows:

1. **Production and startup costs:**
    For each unit `u`, the vectors `u.min_power_cost` and `u.cost_segments`
    are multiplied by a constant `α[u]` sampled from the provided `cost`
    distribution. If `randomize_costs` is false, skips this step.

2. **Load share:**
   For each bus `b` and time `t`, the value `b.load[t]` is multiplied by
   `(β[b] * b.load[t]) / sum(β[b2] * b2.load[t] for b2 in buses)`, where
   `β[b]` is sampled from the provided `load_share` distribution. If
   `randomize_load_share` is false, skips this step.

3. **Peak system load and temporal load profile:**
    Sets the peak load to `ρ * C`, where `ρ` is sampled from `peak_load` and `C`
    is the maximum system capacity, at any time. Also scales the loads of all
    buses, so that `system_load[t+1]` becomes equal to `system_load[t] * γ[t]`,
    where `γ[t]` is sampled from `Normal(load_profile_mu[t], load_profile_sigma[t])`.
    
    The system load for the first time period is set so that the peak load
    matches `ρ * C`. If `load_profile_sigma` and `load_profile_mu` have fewer
    elements than `instance.time`, wraps around. If `randomize_load_profile`
    is false, skips this step.

The default parameters were obtained based on an analysis of publicly available
bid and hourly data from PJM, corresponding to the month of January, 2017. For
more details, see Section 4.2 of the paper.

# References

- **Xavier, Álinson S., Feng Qiu, and Shabbir Ahmed.** *"Learning to solve
  large-scale security-constrained unit commitment problems."* INFORMS Journal
  on Computing 33.2 (2021): 739-756. DOI: 10.1287/ijoc.2020.0976

"""
Base.@kwdef struct Randomization
    cost = Uniform(0.95, 1.05)
    load_profile_mu::Vector{Float64} = [
        1.0,
        0.978,
        0.98,
        1.004,
        1.02,
        1.078,
        1.132,
        1.018,
        0.999,
        1.006,
        0.999,
        0.987,
        0.975,
        0.984,
        0.995,
        1.005,
        1.045,
        1.106,
        0.981,
        0.981,
        0.978,
        0.948,
        0.928,
        0.953,
    ]
    load_profile_sigma::Vector{Float64} = [
        0.0,
        0.011,
        0.015,
        0.01,
        0.012,
        0.029,
        0.055,
        0.027,
        0.026,
        0.023,
        0.013,
        0.012,
        0.014,
        0.011,
        0.008,
        0.008,
        0.02,
        0.02,
        0.016,
        0.012,
        0.014,
        0.015,
        0.017,
        0.024,
    ]
    load_share = Uniform(0.90, 1.10)
    peak_load = Uniform(0.6 * 0.925, 0.6 * 1.075)
    randomize_load_profile::Bool = true
    randomize_costs::Bool = true
    randomize_load_share::Bool = true
end

function _randomize_costs(
    rng,
    sc::ucRHScenario,
    distribution,
)::Nothing
    for unit in sc.thermal_units
        α = rand(rng, distribution)
        unit.min_power_cost *= α
        for k in unit.cost_segments
            k.cost *= α
        end
        for s in unit.startup_categories
            s.cost *= α
        end
    end
    for pu in sc.profiled_units
        α = rand(rng, distribution)
        pu.cost *= α
    end
    for su in sc.storage_units
        α = rand(rng, distribution)
        su.charge_cost *= α
        su.discharge_cost *= α
    end
    return
end

function _randomize_load_share(
    rng,
    sc::ucRHScenario,
    distribution,
)::Nothing
    α = rand(rng, distribution, length(sc.buses))
    for t in 1:sc.time
        total = sum(bus.load[t] for bus in sc.buses)
        den =
            sum(bus.load[t] / total * α[i] for (i, bus) in enumerate(sc.buses))
        for (i, bus) in enumerate(sc.buses)
            bus.load[t] *= α[i] / den
        end
    end
    return
end

function _randomize_load_profile(
    rng,
    sc::ucRHScenario,
    params::Randomization,
)::Nothing
    # Generate new system load
    system_load = [1.0]
    for t in 2:sc.time
        idx = (t - 1) % length(params.load_profile_mu) + 1
        gamma = rand(
            rng,
            Normal(params.load_profile_mu[idx], params.load_profile_sigma[idx]),
        )
        push!(system_load, system_load[t-1] * gamma)
    end
    capacity = sum(maximum(u.max_power) for u in sc.thermal_units)
    peak_load = rand(rng, params.peak_load) * capacity
    system_load = system_load ./ maximum(system_load) .* peak_load

    # Scale bus loads to match the new system load
    prev_system_load = sum(b.load for b in sc.buses)
    for b in sc.buses
        for t in 1:sc.time
            b.load[t] *= system_load[t] / prev_system_load[t]
        end
    end

    return
end

end

"""
    function randomize!(
        instance::ucRH.ucRHInstance,
        method::XavQiuAhm2021.Randomization,
        rng = MersenneTwister(),
    )::Nothing

Randomize costs and loads based on the method described in XavQiuAhm2021.
"""
function randomize!(
    instance::ucRH.ucRHInstance,
    method::XavQiuAhm2021.Randomization;
    rng = MersenneTwister(),
)::Nothing
    for sc in instance.scenarios
        randomize!(sc, method; rng)
    end
    return
end

function randomize!(
    sc::ucRH.ucRHScenario,
    method::XavQiuAhm2021.Randomization;
    rng = MersenneTwister(),
)::Nothing
    if method.randomize_costs
        XavQiuAhm2021._randomize_costs(rng, sc, method.cost)
    end
    if method.randomize_load_share
        XavQiuAhm2021._randomize_load_share(rng, sc, method.load_share)
    end
    if method.randomize_load_profile
        XavQiuAhm2021._randomize_load_profile(rng, sc, method)
    end
    return
end

"""
    function randomize!(
        instance::ucRHInstance;
        method = ucRH.XavQiuAhm2021.Randomization();
        rng = MersenneTwister(),
    )::Nothing

Randomizes instance parameters according to the provided randomization method.

# Example

```julia
instance = ucRH.read_benchmark("matpower/case118/2017-02-01")
ucRH.randomize!(instance)
model = ucRH.build_model(; instance)
```

"""
function randomize!(
    instance::ucRH.ucRHInstance;
    method = XavQiuAhm2021.Randomization(),
    rng = MersenneTwister(),
)::Nothing
    randomize!(instance, method; rng)
    return
end

export randomize!
