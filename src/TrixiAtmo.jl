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
using Printf: @sprintf
using Static: True, False
using StrideArrays: PtrArray
using StaticArrayInterface: static_size
using LinearAlgebra: norm, dot, det
using Reexport: @reexport
using LoopVectorization: @turbo
using Infiltrator

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
export P4estMeshCubedSphere2D, MetricTermsCrossProduct, MetricTermsInvariantCurl,
       MetricTermsExactSpherical, MetricTermsExactCartesian
export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION,
       EARTH_ROTATION_RATE, SECONDS_PER_DAY
export global2contravariant, contravariant2global, spherical2cartesian,
       transform_initial_condition
export initial_condition_gaussian

export examples_dir
end # module TrixiAtmo
