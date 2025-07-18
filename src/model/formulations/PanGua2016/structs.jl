# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

"""
Formulation described in:

    Pan, K., & Guan, Y. (2016). Strong formulations for multistage stochastic
    self-scheduling unit commitment. Operations Research, 64(6), 1482-1498.
    DOI: https://doi.org/10.1287/opre.2016.1520
"""
module PanGua2016

import ..RampingFormulation

struct Ramping <: RampingFormulation end

end
