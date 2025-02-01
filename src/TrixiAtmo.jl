"""
    🌍 TrixiAtmo 🌍

**TrixiAtmo.jl** is a simulation package for atmospheric models based on
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

See also: [trixi-framework/TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl)
"""
module TrixiAtmo

using Trixi
# Import additional symbols that are not exported by Trixi.jl
# using Trixi:
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False
using LinearAlgebra: norm

foo() = true
bar() = false
baz() = Trixi.examples_dir()

include("equations/equations.jl")
export CompressibleMoistEulerEquations2D, CompressibleRainyEulerEquations2D,
       CompressibleRainyEulerPotentialTemperatureEquations2D,
       CompressibleMoistEulerPotentialTemperatureEquations2D,
       CompressibleRainyEulerEquationsExplicit2D

include("callbacks_stage/callbacks_stage.jl")
export NonlinearSolveDG

end # module TrixiAtmo
