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

When using this equation type, functions which are evaluated pointwise, such as fluxes, 
source terms, and initial conditions take in the extra arguments `cache`, `node`, and 
`element`, corresponding to the `cache` field of a `SemidiscretizationHyperbolic`, the node 
index (for tensor-product elements, this should be a tuple of length `NDIMS`), and the 
element index, respectively. To convert an initial condition given in terms of global 
spherical velocity or momentum components to one given in terms of local contravariant 
components, see [`spherical2contravariant`](@ref).
"""
abstract type AbstractCovariantEquations{NDIMS,
                                         NDIMS_AMBIENT,
                                         NVARS} <: AbstractEquations{NDIMS, NVARS} end

@inline nauxvars(equations::AbstractCovariantEquations{NDIMS}) where {NDIMS} = NDIMS^2 +
                                                                               1
@inline nauxvars(equations::AbstractEquations) = 0

@doc raw"""
    spherical2contravariant(initial_condition, ::AbstractCovariantEquations)

Takes in a function with the signature `initial_condition(x, t)` which returns an initial 
condition given in terms of zonal and meridional velocity or momentum components, and 
returns another function with the signature  
`initial_condition_transformed(x, t, equations, cache, node, element)` which returns
the same initial condition with the velocity or momentum vector given in terms of 
contravariant components.
"""
function spherical2contravariant(initial_condition, ::AbstractCovariantEquations)
    function initial_condition_transformed(x, t, aux_vars, equations)
        return spherical2contravariant(initial_condition(x, t), aux_vars, equations)
    end
    return initial_condition_transformed
end
# Numerical flux plus dissipation which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function (numflux::Trixi.FluxPlusDissipation)(u_ll, u_rr,
                                                      aux_vars_ll, aux_vars_rr,
                                                      orientation_or_normal_direction,
                                                      equations::AbstractCovariantEquations)

    # The flux and dissipation need to be defined for the specific equation set
    flux = numflux.numerical_flux(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                  orientation_or_normal_direction, equations)
    diss = numflux.dissipation(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               orientation_or_normal_direction, equations)
    return flux + diss
end

# Central flux which passes node/element indices and cache. 
# We assume that u_ll and u_rr have been transformed into the same local coordinate system.
@inline function Trixi.flux_central(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                    orientation_or_normal_direction,
                                    equations::AbstractCovariantEquations)
    flux_ll = Trixi.flux(u_ll, aux_vars_ll, orientation_or_normal_direction, equations)
    flux_rr = Trixi.flux(u_rr, aux_vars_rr, orientation_or_normal_direction, equations)
    return 0.5f0 * (flux_ll + flux_rr)
end

# Extract geometric information from auxiliary variables
@inline volume_element(aux_vars, ::AbstractCovariantEquations{2}) = aux_vars[1]
@inline basis_covariant(aux_vars, ::AbstractCovariantEquations{2}) = SMatrix{2, 2}(aux_vars[2:5])

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end

include("reference_data.jl")
include("covariant_advection.jl")
include("compressible_moist_euler_2d_lucas.jl")
include("shallow_water_3d.jl")
end # @muladd
