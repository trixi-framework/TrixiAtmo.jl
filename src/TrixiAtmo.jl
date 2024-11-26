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
using StrideArrays: PtrArray
using StaticArrayInterface: static_size
using LinearAlgebra: norm, dot
using Reexport: @reexport
using LoopVectorization: @turbo

@reexport using StaticArrays: SVector, SMatrix

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("meshes/meshes.jl")
include("semidiscretization/semidiscretization.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization_hyperbolic_2d_manifold_in_3d.jl")
include("callbacks_step/callbacks_step.jl")

export CompressibleMoistEulerEquations2D, ShallowWaterEquations3D,
       CovariantLinearAdvectionEquation2D

export flux_chandrashekar, flux_LMARS

export velocity, waterheight, pressure, energy_total, energy_kinetic, energy_internal,
       lake_at_rest_error, source_terms_lagrange_multiplier,
       clean_solution_lagrange_multiplier!
export P4estCubedSphere2D, P4estMeshQuadIcosahedron2D, MetricTermsCrossProduct,
       MetricTermsInvariantCurl
export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION,
       EARTH_ROTATION_RATE, SECONDS_PER_DAY
export spherical2contravariant, contravariant2spherical, spherical2cartesian,
       transform_to_cartesian, transform_to_contravariant
export initial_condition_gaussian

export examples_dir
end # module TrixiAtmo
