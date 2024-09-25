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
@reexport using StaticArrays: SVector

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("meshes/meshes.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization_hyperbolic_2d_manifold_in_3d.jl")

export CompressibleMoistEulerEquations2D
export CompressibleEulerPotentialTemperatureEquations2D

export flux_chandrashekar, flux_theta, flux_theta_rhoAM, flux_theta_physentropy,
       flux_theta_physentropy2

export examples_dir

end # module TrixiAtmo
