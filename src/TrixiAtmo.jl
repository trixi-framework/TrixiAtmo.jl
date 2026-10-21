"""
    TrixiAtmo

**TrixiAtmo.jl** is a simulation package for atmospheric models based on
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

See also: [trixi-framework/TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl)
"""
module TrixiAtmo

using Accessors: setproperties
using Reexport: @reexport
using Trixi
using MuladdMacro: @muladd
using Printf: @sprintf
using Static: True, False
using StrideArrays: PtrArray
using StaticArrayInterface: static_size
using LinearAlgebra: Diagonal, I, UpperTriangular, cross, norm, dot, det, diagm
using Reexport: @reexport
using LoopVectorization: @turbo
using QuadGK: quadgk
using ForwardDiff: derivative
using HDF5: HDF5, h5open, attributes, create_dataset, datatype, dataspace

@reexport using StaticArrays: SVector, SMatrix
@reexport import Trixi: waterheight, varnames, cons2cons, cons2prim,
                        prim2cons, cons2entropy, entropy2cons, velocity,
                        max_abs_speeds, max_abs_speed_naive, max_abs_speed,
                        have_nonconservative_terms, boundary_condition_slip_wall,
                        energy_kinetic, energy_internal, energy_total, entropy, pressure,
                        flux, flux_ec, flux_chandrashekar, flux_wintermeyer_etal,
                        flux_fjordholm_etal, flux_nonconservative_wintermeyer_etal,
                        flux_nonconservative_fjordholm_etal, FluxLMARS, flux_shima_etal,
                        flux_ranocha, flux_kennedy_gruber, trixi_backend

using Trixi: ln_mean, stolarsky_mean, inv_ln_mean

# Time integration
using Trixi: AbstractTimeIntegrationAlgorithm, TimeIntegratorSolution,
             @threaded, @trixi_timeit, @trixi_timeit_ext, timer
using SciMLBase: SciMLBase, CallbackSet, ODEProblem
import SciMLBase: init, step!, add_tstop!, get_tmp_cache, terminate!, get_du,
                  get_proposed_dt, set_proposed_dt!, has_tstop, first_tstop
using DiffEqBase: DiffEqBase, get_tstops, get_tstops_array
# To keep backwards compatibility with SciMLBase v2 as in Trixi.jl, see
# https://github.com/trixi-framework/Trixi.jl/pull/2918#issuecomment-4233720339
@static if isdefined(SciMLBase, :derivative_discontinuity!)
    import SciMLBase: derivative_discontinuity!
else
    const derivative_discontinuity! = SciMLBase.u_modified!
end
using DataStructures: BinaryHeap, FasterForward
using KernelAbstractions: @kernel, @index, Backend

# DGMulti Solvers
using StartUpDG: MeshData, RefElemData

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("meshes/meshes.jl")
include("semidiscretization/semidiscretization.jl")
include("solvers/solvers.jl")
include("callbacks_step/callbacks_step.jl")
include("callbacks_stage/callbacks_stage.jl")
include("time_integration/time_integration.jl")

export CompressibleMoistEulerEquations2D,
       CompressibleRainyEulerEquations2D,
       CovariantLinearAdvectionEquation2D, CovariantShallowWaterEquations2D,
       CovariantEulerEnergyEquations2D,
       ShallowWaterEquations3D, SplitCovariantShallowWaterEquations2D,
       CompressibleEulerPotentialTemperatureEquations1D,
       CompressibleEulerPotentialTemperatureEquations2D,
       CompressibleEulerPotentialTemperatureEquations3D,
       CompressibleEulerPotentialTemperatureEquationsWithGravity1D,
       CompressibleEulerPotentialTemperatureEquationsWithGravity2D,
       CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
       CompressibleEulerEnergyEquationsWithGravity2D,
       CompressibleEulerEnergyEquationsWithGravity3D,
       CompressibleEulerInternalEnergyEquationsWithGravity2D

export GlobalCartesianCoordinates, GlobalSphericalCoordinates

export NonlinearSolveDG

export SSPRW32, ROS2, ARK2, ARS232, ARS343

export flux_chandrashekar, FluxLMARS

export flux_nonconservative_zeros, flux_nonconservative_ec,
       flux_nonconservative_surface_simplified, source_terms_geometric_coriolis,
       source_terms_coriolis, source_terms_coriolis_lagrange_multiplier,
       source_terms_gravity, flux_tec, flux_etec, flux_nonconservative_souza_etal,
       flux_nonconservative_artiano_etal,
       flux_nonconservative_waruszewski_etal, flux_zero,
       flux_ec_rain, flux_LMARS, flux_nonconservative_es, flux_conservative_es,
       flux_conservative_etec, flux_nonconservative_etec

export flux_lmars_combined, flux_kennedy_gruber_souza_etal, flux_lmars_es

export source_terms_lagrange_multiplier, clean_solution_lagrange_multiplier!

export cons2prim_and_vorticity, contravariant2global, contravariant_cons2global_prim

export P4estMeshCubedSphere2D, P4estMeshQuadIcosahedron2D, DGMultiMeshTriIcosahedron2D,
       MetricTermsCrossProduct, MetricTermsInvariantCurl, MetricTermsCovariant,
       SphericalManifold, FlatManifold,
       ChristoffelSymbolsAutodiff, ChristoffelSymbolsCollocationDerivative

export EARTH_RADIUS, EARTH_GRAVITATIONAL_ACCELERATION,
       EARTH_ROTATION_RATE, SECONDS_PER_DAY,
       SMALL_EARTH_SCALING, SMALL_EARTH_RADIUS, SMALL_EARTH_ROTATION_RATE

export transform_initial_condition

export initial_condition_gaussian, initial_condition_geostrophic_balance,
       initial_condition_rossby_haurwitz, initial_condition_isolated_mountain,
       initial_condition_unsteady_solid_body_rotation,
       initial_condition_barotropic_instability

export initial_condition_adiabatic, initial_condition_gravity_waves,
       initial_condition_baroclinic_instability, initial_condition_isothermal,
       initial_condition_vortex_shedding

export boundary_condition_horizontal_slip_wall

export source_terms_coriolis_small_earth

export MountainWaveSetup, SchärMountain, WitchOfAgnesi, terrain_following_faces

export cartesian_to_sphere, rayleigh_damping_coefficient

export bottom_topography_isolated_mountain, bottom_topography_unsteady_solid_body_rotation

export AtmosphereLayers, AtmosphereLayersRainyBubble

export P4estMeshCubedSphereTopography, GalChen, Sleve

export examples_dir

end # module TrixiAtmo
