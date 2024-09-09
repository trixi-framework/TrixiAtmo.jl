using Trixi: AbstractEquations

# Physical constants
const EARTH_RADIUS = 6.37122f6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616f0  # m/s²
const EARTH_ROTATION_RATE = 7.292f-5  # rad/s
const SECONDS_PER_DAY = 8.64f4

# Abstract type used to dispatch specialized solvers for the covariant form
abstract type AbstractCovariantEquations{NDIMS, NVARS} <: AbstractEquations{NDIMS, NVARS} end

# Numerical flux plus dissipation which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function (numflux::Trixi.FluxPlusDissipation)(u_ll, u_rr,
                                                      orientation_or_normal_direction,
                                                      equations::AbstractCovariantEquations,
                                                      i, j, element, cache)

    # The flux and dissipation need to be defined for the specific equation set
    flux = numflux.numerical_flux(u_ll, u_rr, orientation_or_normal_direction, equations,
                                  i, j, i, j, element, cache)
    diss = numflux.dissipation(u_ll, u_rr, orientation_or_normal_direction, equations,
                               i, j, element, cache)
    return flux + diss
end

# Central flux which passes node/element indices and cache. 
@inline function Trixi.flux_central(u_ll, u_rr,
                                    orientation_or_normal_direction,
                                    equations::AbstractCovariantEquations,
                                    i_ll, j_ll, i_rr, j_rr, element, cache)
    flux_ll = Trixi.flux(u_ll, orientation_or_normal_direction, equations,
                         i_ll, j_ll, element, cache)
    flux_rr = Trixi.flux(u_rr, orientation_or_normal_direction, equations,
                         i_rr, j_rr, element, cache)

    return 0.5f0 * (flux_ll + flux_rr)
end

include("covariant_advection.jl")
abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
include("compressible_moist_euler_2d_lucas.jl")
