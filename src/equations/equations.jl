using Trixi: AbstractEquations

# Physical constants
const EARTH_RADIUS = 6.37122f6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616f0  # m/s²
const EARTH_ROTATION_RATE = 7.292f-5  # rad/s
const SECONDS_PER_DAY = 8.64f4

# Abstract type used to dispatch specialized solvers for the covariant form
abstract type AbstractCovariantEquations{NDIMS,
                                         NDIMS_AMBIENT,
                                         NVARS} <: AbstractEquations{NDIMS, NVARS} end

# Numerical flux plus dissipation which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function (numflux::Trixi.FluxPlusDissipation)(u_ll, u_rr,
                                                      orientation_or_normal_direction,
                                                      equations::AbstractCovariantEquations,
                                                      elements, i, j, element)

    # The flux and dissipation need to be defined for the specific equation set
    flux = numflux.numerical_flux(u_ll, u_rr, orientation_or_normal_direction, equations,
                                  elements, i, j, i, j, element)
    diss = numflux.dissipation(u_ll, u_rr, orientation_or_normal_direction, equations,
                               elements, i, j, element)
    return flux + diss
end

# Central flux which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function Trixi.flux_central(u_ll, u_rr,
                                    orientation_or_normal_direction,
                                    equations::AbstractCovariantEquations,
                                    elements, i_ll, j_ll, i_rr, j_rr, element)
    flux_ll = Trixi.flux(u_ll, orientation_or_normal_direction, equations,
                         elements, i_ll, j_ll, element)
    flux_rr = Trixi.flux(u_rr, orientation_or_normal_direction, equations,
                         elements, i_rr, j_rr, element)

    return 0.5f0 * (flux_ll + flux_rr)
end

include("covariant_advection.jl")
include("covariant_shallow_water.jl")
abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
include("compressible_moist_euler_2d_lucas.jl")
