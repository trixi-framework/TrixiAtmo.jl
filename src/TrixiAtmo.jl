"""
    TrixiAtmo

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
using QuadGK: quadgk
using ForwardDiff: derivative
using HDF5: HDF5, h5open, attributes, create_dataset, datatype, dataspace

@reexport using StaticArrays: SVector, SMatrix
@reexport import Trixi: waterheight

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("callbacks_stage/callbacks_stage.jl")
include("meshes/meshes.jl")
include("semidiscretization/semidiscretization.jl")
include("solvers/solvers.jl")
include("callbacks_step/callbacks_step.jl")

export CompressibleMoistEulerEquations2D, ShallowWaterEquations3D,
       CovariantLinearAdvectionEquation2D, CovariantShallowWaterEquations2D,
       SplitCovariantShallowWaterEquations2D, CompressibleRainyEulerEquations2D,
       CompressibleMoistEulerPotentialTemperatureEquations2D,
       CompressibleRainyEulerExplicitEquations2D

export GlobalCartesianCoordinates, GlobalSphericalCoordinates

export NonlinearSolveDG

export flux_chandrashekar, FluxLMARS

export flux_nonconservative_zeros, flux_nonconservative_ec,
       flux_nonconservative_surface_simplified, source_terms_geometric_coriolis,
       source_terms_coriolis, source_terms_coriolis_lagrange_multiplier

export source_terms_lagrange_multiplier, clean_solution_lagrange_multiplier!

export cons2prim_and_vorticity, contravariant2global

export P4estMeshCubedSphere2D, P4estMeshQuadIcosahedron2D, MetricTermsCrossProduct,
       MetricTermsInvariantCurl, MetricTermsCovariantSphere, ChristoffelSymbolsAutodiff,
       ChristoffelSymbolsCollocationDerivative

export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION,
       EARTH_ROTATION_RATE, SECONDS_PER_DAY

export transform_initial_condition

export initial_condition_gaussian, initial_condition_geostrophic_balance,
       initial_condition_rossby_haurwitz, initial_condition_isolated_mountain,
       initial_condition_unsteady_solid_body_rotation,
       initial_condition_barotropic_instability

export bottom_topography_isolated_mountain, bottom_topography_unsteady_solid_body_rotation

export examples_dir
end # module TrixiAtmo
