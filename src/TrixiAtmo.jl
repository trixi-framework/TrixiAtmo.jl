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
using Static: True, False
using StaticArrayInterface: static_size
using StrideArrays: StrideArray, StaticInt, PtrArray
using LinearAlgebra: norm, dot
using DiffEqCallbacks: PeriodicCallback, PeriodicCallbackAffect

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("meshes/meshes.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization_hyperbolic_2d_manifold_in_3d.jl")

export CompressibleMoistEulerEquations2D

export flux_chandrashekar, flux_LMARS

export examples_dir
end # module TrixiAtmo
