"""
    🌍 TrixiAtmo 🌍

**TrixiAtmo.jl** is a simulation package for atmospheric models based on
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

See also: [trixi-framework/TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl)
"""
module TrixiAtmo

using Trixi
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False
using LinearAlgebra: norm
using Reexport: @reexport
@reexport using StaticArrays: SVector

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")

export CompressibleMoistEulerEquations2D

export flux_chandrashekar, flux_LMARS

export examples_dir

end # module TrixiAtmo
