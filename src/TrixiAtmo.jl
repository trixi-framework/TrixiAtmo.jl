"""
    🌍 TrixiAtmo 🌍

**TrixiAtmo.jl** is a simulation package for atmospheric models based on
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

See also: [trixi-framework/TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl)
"""
module TrixiAtmo

using Trixi
using MuladdMacro: @muladd
using Static: True, False
using StrideArrays: PtrArray
using StaticArrayInterface: static_size
using LinearAlgebra: norm
using Reexport: @reexport
using LoopVectorization: @turbo
@reexport using StaticArrays: SVector

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("meshes/meshes.jl")
include("semidiscretization/semidiscretization.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization_hyperbolic_2d_manifold_in_3d.jl")
include("callbacks_step/stepsize_dg2d.jl")

export CompressibleMoistEulerEquations2D, ShallowWaterEquations3D

export flux_chandrashekar, flux_LMARS

export velocity, waterheight, pressure, energy_total, energy_kinetic, energy_internal,
       lake_at_rest_error

export examples_dir

end # module TrixiAtmo
