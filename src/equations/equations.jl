using Trixi: AbstractEquations

# Physical constants
const EARTH_RADIUS = 6.37122f6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616f0  # m/s²
const EARTH_ROTATION_RATE = 7.292f-5  # rad/s
const SECONDS_PER_DAY = 8.64f4

# Abstract type used to dispatch specialized solvers for the covariant form of a partial
# differential equation on a two-dimensional manifold in 3D space
abstract type AbstractCovariantEquations2D{NVARS} <: AbstractEquations{2, NVARS} end
include("covariant_advection.jl")

# Numerical flux which passes node/element indices and cache. We assume that u_ll and u_rr have 
# been transformed into the same local coordinate system.
@inline function (numflux::Trixi.FluxPlusDissipation)(u_ll, u_rr,
                                                      orientation_or_normal_direction,
                                                      equations::AbstractCovariantEquations2D,
                                                      i, j, element, cache)

    # The flux and dissipation need to be defined for the specific equation set
    flux = numflux.numerical_flux(u_ll, u_rr, orientation_or_normal_direction, equations,
                                  i, j, i, j, element, cache)
    diss = numflux.dissipation(u_ll, u_rr, orientation_or_normal_direction, equations,
                               i, j, element, cache)
    return flux + diss
end

@inline function Trixi.flux_central(u_ll, u_rr,
                                    orientation_or_normal_direction,
                                    equations::AbstractCovariantEquations2D,
                                    i_ll, j_ll, i_rr, j_rr, element, cache)
    flux_ll = Trixi.flux(u_ll, orientation_or_normal_direction, equations,
                         i_ll, j_ll, element, cache)
    flux_rr = Trixi.flux(u_rr, orientation_or_normal_direction, equations,
                         i_rr, j_rr, element, cache)

    return 0.5f0 * (flux_ll + flux_rr)
end

# Convert from spherical polar coordinates to contravariant components.
# We assume the first component of u_pol is a conserved scalar
# and the second and third define spherical vector components to be transformed
@inline function pol2con(u_pol::SVector{3}, ::AbstractCovariantEquations2D, i, j, element,
                         cache)
    A11, A21, A12, A22 = view(cache.elements.polar_transform_matrix, :, :, i, j, element)
    inv_detA = 1 / (A11 * A22 - A12 * A21)
    return SVector(u_pol[1], (A22 * u_pol[2] - A12 * u_pol[3]) * inv_detA,
                   (-A21 * u_pol[2] + A11 * u_pol[3]) * inv_detA)
end

# Convert from spherical polar coordinates to contravariant components.
# We assume the first component of u_con is a conserved scalar
# and the second and third define contravariant vector components to be transformed
@inline function con2pol(u_con::SVector{3}, ::AbstractCovariantEquations2D, i, j, element,
                         cache)
    A11, A21, A12, A22 = view(cache.elements.polar_transform_matrix, :, :, i, j, element)
    return SVector(u_con[1], A11 * u_con[2] + A12 * u_con[3],
                   A21 * u_con[2] + A22 * u_con[3])
end

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
include("compressible_moist_euler_2d_lucas.jl")
