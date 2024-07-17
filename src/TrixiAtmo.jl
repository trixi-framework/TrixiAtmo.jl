"""
    🌍 TrixiAtmo 🌍

**TrixiAtmo.jl** is a simulation package for atmospheric models based on
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

See also: [trixi-framework/TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl)
"""
module TrixiAtmo

using Reexport: @reexport

using Trixi
using MuladdMacro: @muladd
@reexport using StaticArrays: SVector, SMatrix
using Static: True, False
using LinearAlgebra: norm, dot
using DiffEqCallbacks: PeriodicCallback, PeriodicCallbackAffect

export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION, EARTH_ROTATION_RATE, SECONDS_PER_DAY
include("equations/equations.jl")

include("solvers/solvers.jl")

end # module TrixiAtmo
