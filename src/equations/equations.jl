@muladd begin
#! format: noindent

using Trixi: AbstractEquations

@doc raw"""
    AbstractCovariantEquations{NDIMS, 
                               NDIMS_AMBIENT, 
                               GlobalCoordinateSystem,
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
local contravariant components, see [`transform_to_contravariant`](@ref).
"""
abstract type AbstractCovariantEquations{NDIMS,
                                         NDIMS_AMBIENT,
                                         GlobalCoordinateSystem,
                                         NVARS} <: AbstractEquations{NDIMS, NVARS} end

# Coordinate systems
struct GlobalCartesianCoordinates end
struct GlobalSphericalCoordinates end

"""
    have_aux_node_vars(equations)

Trait function determining whether `equations` requires the use of auxiliary variables.
Classical conservation laws such as the [`CompressibleEulerEquations2D`](@ref) do not 
require auxiliary variables. The return value will be `True()` or `False()` to allow 
dispatching on the return type.
"""
@inline have_aux_node_vars(::AbstractEquations) = False()

# cons2cons method which takes in unused aux_vars variable
@inline Trixi.cons2cons(u, aux_vars, equations) = u

# If no auxiliary variables are passed into the conversion to spherical coordinates, do not 
# do any conversion.
@inline contravariant2global(u, equations) = u

# For the covariant form, the auxiliary variables are the the NDIMS*NDIMS_AMBIENT entries 
# of the covariant basis matrix
@inline have_aux_node_vars(::AbstractCovariantEquations) = True()
@inline n_aux_node_vars(::AbstractCovariantEquations{NDIMS, NDIMS_AMBIENT}) where {NDIMS,
NDIMS_AMBIENT} = 2 * NDIMS * NDIMS_AMBIENT + 1

# Return auxiliary variable names for 2D covariant form
@inline function auxvarnames(::AbstractCovariantEquations{2})
    return ("basis_covariant[1,1]",
            "basis_covariant[2,1]",
            "basis_covariant[3,1]",
            "basis_covariant[1,2]",
            "basis_covariant[2,2]",
            "basis_covariant[3,2]",
            "basis_contravariant[1,1]",
            "basis_contravariant[2,1]",
            "basis_contravariant[3,1]",
            "basis_contravariant[1,2]",
            "basis_contravariant[2,2]",
            "basis_contravariant[3,2]",
            "area_element")
end

# Extract the covariant basis vectors a_i from the auxiliary variables as a matrix A, 
# where a_i = A[:, i] denotes the covariant tangent basis in a spherical/ellipsoidal 
# coordinate system.
@inline function basis_covariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{3, 2}(aux_vars[1], aux_vars[2], aux_vars[3],
                         aux_vars[4], aux_vars[5], aux_vars[6])
end

@inline function basis_contravariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{2, 3}(aux_vars[7], aux_vars[8],
                         aux_vars[9], aux_vars[10],
                         aux_vars[11], aux_vars[12])
end

@inline function area_element(aux_vars, ::AbstractCovariantEquations{2})
    return aux_vars[13]
end

# Transform zonal and meridional velocity/momentum components to Cartesian components
function spherical2cartesian(vlon, vlat, x)
    # Co-latitude
    colat = acos(x[3] / sqrt(x[1]^2 + x[2]^2 + x[3]^2))

    # Longitude
    if sign(x[2]) == 0.0
        signy = 1.0
    else
        signy = sign(x[2])
    end
    r_xy = sqrt(x[1]^2 + x[2]^2)
    if r_xy == 0.0
        lon = pi / 2
    else
        lon = signy * acos(x[1] / r_xy)
    end

    v1 = -cos(colat) * cos(lon) * vlat - sin(lon) * vlon
    v2 = -cos(colat) * sin(lon) * vlat + cos(lon) * vlon
    v3 = sin(colat) * vlat

    return SVector(v1, v2, v3)
end

@doc raw"""
    transform_initial_condition(initial_condition, equations)

Takes in a function with the signature `initial_condition(x, t)` which returns an initial 
condition given in terms of global velocity or momentum components, and returns another
function with the signature  `initial_condition_transformed(x, t, aux_vars, equations)` 
which returns the same initial condition with the velocity or momentum vector given in
terms of contravariant components.
"""
function transform_initial_condition(initial_condition, ::AbstractCovariantEquations)
    function initial_condition_transformed(x, t, aux_vars, equations)
        return global2contravariant(initial_condition(x, t), aux_vars, equations)
    end
    return initial_condition_transformed
end

function transform_initial_condition(initial_condition, ::AbstractEquations)
    function initial_condition_transformed(x, t, equations)
        return Trixi.prim2cons(initial_condition(x, t), equations)
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
