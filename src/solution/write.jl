# ucRH.jl: Optimization Package for Security-Constrained Unit Commitment
# Copyright (C) 2020, UChicago Argonne, LLC. All rights reserved.
# Released under the modified BSD license. See COPYING.md for more details.

"""
    write(filename::AbstractString, solution::AbstractDict)::Nothing

Write the given solution to a JSON file.

# Example

```julia
solution = ucRH.solution(model)
ucRH.write("/tmp/output.json", solution)
```
"""
function write(filename::AbstractString, solution::AbstractDict)::Nothing
    open(filename, "w") do file
        return JSON.print(file, solution, 2)
    end
    return
end
