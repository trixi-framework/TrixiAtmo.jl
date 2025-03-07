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
using LinearAlgebra: cross, norm, dot, det
using Reexport: @reexport
using LoopVectorization: @turbo
using ForwardDiff: derivative
using HDF5: HDF5, h5open, attributes, create_dataset, datatype, dataspace

@reexport using StaticArrays: SVector, SMatrix

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("meshes/meshes.jl")
include("semidiscretization/semidiscretization.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization_hyperbolic_2d_manifold_in_3d.jl")
include("callbacks_step/callbacks_step.jl")

export CompressibleMoistEulerEquations2D, ShallowWaterEquations3D,
       CovariantLinearAdvectionEquation2D, CovariantShallowWaterEquations2D,
       SplitCovariantShallowWaterEquations2D

export GlobalCartesianCoordinates, GlobalSphericalCoordinates

export flux_chandrashekar, FluxLMARS

export flux_nonconservative_zeros, flux_nonconservative_ec,
       source_terms_geometric_coriolis

export velocity, waterheight, pressure, energy_total, energy_kinetic, energy_internal,
       lake_at_rest_error, source_terms_lagrange_multiplier,
       clean_solution_lagrange_multiplier!

export cons2prim_and_vorticity

export P4estMeshCubedSphere2D, P4estMeshQuadIcosahedron2D, MetricTermsCrossProduct,
       MetricTermsInvariantCurl

export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION,
       EARTH_ROTATION_RATE, SECONDS_PER_DAY

export global2contravariant, contravariant2global, spherical2cartesian, cartesian2spherical,
       transform_initial_condition

export initial_condition_gaussian, initial_condition_geostrophic_balance,
       initial_condition_rossby_haurwitz

export examples_dir
end # module TrixiAtmo
