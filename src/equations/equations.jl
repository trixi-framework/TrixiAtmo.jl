@muladd begin
#! format: noindent

using Trixi: AbstractEquations

@doc raw"""
    AbstractCovariantEquations{NDIMS, 
                               NDIMS_AMBIENT, 
                               NVARS} <: AbstractEquations{NDIMS, NVARS} 

Abstract type used to dispatch on systems of equations in covariant form, in which fluxes
and prognostic variables are stored and computed in terms of their contravariant components 
defining their expansions in terms of the local covariant tangent basis. The type parameter
`NDIMS` denotes the dimension of the manifold on which the equations are solved, while
`NDIMS_AMBIENT` is the dimension of the ambient space in which such a manifold is embedded. 
Some references on discontinuous Galerkin methods in covariant flux form are listed below:

- M. Baldauf (2020). Discontinuous Galerkin solver for the shallow-water equations in
  covariant form on the sphere and the ellipsoid. Journal of Computational Physics 
  410:109384. [DOI: 10.1016/j.jcp.2020.109384](https://doi.org/10.1016/j.jcp.2020.109384) 
- M. Baldauf (2021). A horizontally explicit, vertically implicit (HEVI) discontinuous
  Galerkin scheme for the 2-dimensional Euler and Navier-Stokes equations using 
  terrain-following coordinates. Journal of Computational Physics 446:110635. [DOI: 10.1016/
  j.jcp.2021.110635](https://doi.org/10.1016/j.jcp.2021.110635)
- L. Bao, R. D. Nair, and H. M. Tufo (2014). A mass and momentum flux-form high-order
  discontinuous Galerkin shallow water model on the cubed-sphere. A mass and momentum 
  flux-form high-order discontinuous Galerkin shallow water model on the cubed-sphere. 
  Journal of Computational Physics 271:224-243. 
  [DOI: 10.1016/j.jcp.2013.11.033](https://doi.org/10.1016/j.jcp.2013.11.033)

!!! note 
    Components of vector-valued fields should be prescibed within the global coordinate 
    system (i.e. zonal and meridional components in the case of a spherical shell). 
    When dispatched on this type, the function `Trixi.compute_coefficients!` will 
    internally use the `covariant_basis` field of the container type 
    [`P4estElementContainerCovariant`](@ref) to obtain the local contravariant components
    used in the solver. 
"""
abstract type AbstractCovariantEquations{NDIMS,
                                         NDIMS_AMBIENT,
                                         NVARS} <: AbstractEquations{NDIMS, NVARS} end

# Numerical flux plus dissipation which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function (numflux::Trixi.FluxPlusDissipation)(u_ll, u_rr,
                                                      orientation_or_normal_direction,
                                                      equations::AbstractCovariantEquations{2},
                                                      cache, node_ll, node_rr, element)

    # The flux and dissipation need to be defined for the specific equation set
    flux = numflux.numerical_flux(u_ll, u_rr, orientation_or_normal_direction,
                                  equations, cache, node_ll, node_rr, element)
    diss = numflux.dissipation(u_ll, u_rr, orientation_or_normal_direction, equations,
                               cache, node_ll, node_rr, element)
    return flux + diss
end

# Central flux which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function Trixi.flux_central(u_ll, u_rr,
                                    orientation_or_normal_direction,
                                    equations::AbstractCovariantEquations{2},
                                    cache, node_ll, node_rr, element)
    flux_ll = Trixi.flux(u_ll, orientation_or_normal_direction, equations,
                         cache, node_ll, element)
    flux_rr = Trixi.flux(u_rr, orientation_or_normal_direction, equations,
                         cache, node_rr, element)

    return 0.5f0 * (flux_ll + flux_rr)
end

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end

include("reference_data.jl")
include("covariant_advection.jl")
include("compressible_moist_euler_2d_lucas.jl")
include("shallow_water_3d.jl")
end # @muladd
