module TrixiAtmo

using Reexport: @reexport

using Trixi
using MuladdMacro: @muladd
@reexport using StaticArrays: SVector, SMatrix
using Static: True, False
using LinearAlgebra: norm, dot
using DiffEqCallbacks: PeriodicCallback, PeriodicCallbackAffect

include("meshes/meshes.jl")

export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION, EARTH_ROTATION_RATE, SECONDS_PER_DAY
include("equations/equations.jl")

end # module TrixiAtmo