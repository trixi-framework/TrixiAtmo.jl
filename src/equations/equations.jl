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
source terms, and initial conditions take in the extra argument `aux_vars`, which contains 
the geometric information needed for the covariant form. To convert an initial condition 
given in terms of global spherical velocity or momentum components to one given in terms of 
local contravariant components, see [`spherical2contravariant`](@ref).
"""
abstract type AbstractCovariantEquations{NDIMS,
                                         NDIMS_AMBIENT,
                                         NVARS} <: AbstractEquations{NDIMS, NVARS} end

"""
    have_aux_node_vars(equations)

Trait function determining whether `equations` requires the use of auxiliary variables.
Classical conservation laws such as the [`CompressibleEulerEquations2D`](@ref) do not 
require auxiliary variables. The return value will be `True()` or `False()` to allow 
dispatching on the return type.
"""
@inline have_aux_node_vars(::AbstractEquations) = False()

# For the covariant form, the auxiliary variables are the the NDIMS^2 entries of the 
# covariant basis matrix
@inline have_aux_node_vars(::AbstractCovariantEquations) = True()
@inline n_aux_node_vars(::AbstractCovariantEquations{NDIMS}) where {NDIMS} = NDIMS^2

# Return auxiliary variable names for 2D covariant form
@inline function auxvarnames(::AbstractCovariantEquations{2})
    return ("basis_covariant[1,1]", "basis_covariant[2,1]",
            "basis_covariant[1,2]", "basis_covariant[2,2]")
end

# Extract the covariant basis vectors a_i from the auxiliary variables as a matrix A, 
# where a_i = A[:, i] denotes the covariant tangent basis in a spherical/ellipsoidal 
# coordinate system.
@inline function basis_covariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{2, 2}(aux_vars[1], aux_vars[2], aux_vars[3], aux_vars[4])
end
@inline function area_element(aux_vars, ::AbstractCovariantEquations{2})
    return abs(aux_vars[1] * aux_vars[4] - aux_vars[2] * aux_vars[3])
end

@doc raw"""
    spherical2contravariant(initial_condition, ::AbstractCovariantEquations)

Takes in a function with the signature `initial_condition(x, t)` which returns an initial 
condition given in terms of zonal and meridional velocity or momentum components, and 
returns another function with the signature  `initial_condition_transformed(x, t, aux_vars, 
equations)` which returns the same initial condition with the velocity or momentum vector given in terms of contravariant components.
"""
function spherical2contravariant(initial_condition, ::AbstractCovariantEquations)
    function initial_condition_transformed(x, t, aux_vars, equations)
        return spherical2contravariant(initial_condition(x, t), aux_vars, equations)
    end
    return initial_condition_transformed
end

# Numerical flux plus dissipation for abstract covariant equations as a function of the 
# state vectors u_ll and u_rr, as well as the auxiliary variable vectors aux_vars_ll and 
# aux_vars_rr, which contain the geometric information. We assume that u_ll and u_rr have 
# been transformed into the same local coordinate system.
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

# Central flux for abstract covariant equations as a function of the state vectors u_ll and 
# u_rr, as well as the auxiliary variable vectors aux_vars_ll and aux_vars_rr, which 
# contain the geometric information. We assume that u_ll and u_rr have been transformed 
# into the same  local coordinate system.
@inline function Trixi.flux_central(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                    orientation_or_normal_direction,
                                    equations::AbstractCovariantEquations)
    flux_ll = Trixi.flux(u_ll, aux_vars_ll, orientation_or_normal_direction, equations)
    flux_rr = Trixi.flux(u_rr, aux_vars_rr, orientation_or_normal_direction, equations)
    return 0.5f0 * (flux_ll + flux_rr)
end

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end

include("reference_data.jl")
include("covariant_advection.jl")
include("compressible_moist_euler_2d_lucas.jl")
include("shallow_water_3d.jl")
end # @muladd
