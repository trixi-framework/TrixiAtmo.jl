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
the geometric information needed for the covariant form. The type parameter 
`GlobalCoordinateSystem` specifies the global coordinate system used to define the 
covariant tangent basis, and may be either [`GlobalCartesianCoordinates`](@ref) or 
[`GlobalSphericalCoordinates`](@ref). The `GlobalCoordinateSystem` type parameter also 
specifies the coordinate system with respect to which the initial condition should be 
prescribed.
"""
abstract type AbstractCovariantEquations{NDIMS,
                                         NDIMS_AMBIENT,
                                         GlobalCoordinateSystem,
                                         NVARS} <: AbstractEquations{NDIMS, NVARS} end

"""
    GlobalCartesianCoordinates()

Struct used for dispatch, specifying that the covariant tangent basis vectors should be 
defined with respect to a global Cartesian coordinate system.
"""
struct GlobalCartesianCoordinates end

"""
    GlobalSphericalCoordinates()

Struct used for dispatch, specifying that the covariant tangent basis vectors should be 
defined with respect to a global spherical coordinate system.
"""
struct GlobalSphericalCoordinates end

"""
    have_aux_node_vars(equations)

Trait function determining whether `equations` requires the use of auxiliary variables.
Classical conservation laws such as the [`CompressibleEulerEquations2D`](@ref) do not 
require auxiliary variables. The return value will be `True()` or `False()` to allow 
dispatching on the return type.
"""
@inline have_aux_node_vars(::AbstractEquations) = False()

@doc raw"""
    transform_initial_condition(initial_condition, equations)

Takes in a function with the signature `initial_condition(x, t, equations)` which returns 
an initial condition given in terms of global Cartesian or zonal/meridional velocity
components, and returns another function `initial_condition_transformed(x, t, equations)` 
or `initial_condition_transformed(x, t, aux_vars, equations)` which returns the same 
initial data, but transformed to the appropriate prognostic variables used internally by 
the solver. For the covariant form, this involves a transformation of the global velocity 
components to contravariant components using [`global2contravariant`](@ref) as well as a 
conversion from primitive to conservative variables. For standard Cartesian formulations, 
this simply involves a conversion from  primitive to conservative variables. The intention 
here is to have a set of test cases (for example, [`initial_condition_gaussian`](@ref)) for 
which the initial condition is prescribed using a standardized set of primitive variables 
in a global coordinate system, and transformed to the specific prognostic variables 
required for a given model.
!!! note 
    When using the covariant formulation, the initial velocity components should be defined 
    in the coordinate system specified by the `GlobalCoordinateSystem` type parameter in
    [`AbstractCovariantEquations`](@ref).
"""
function transform_initial_condition(initial_condition, ::AbstractCovariantEquations)
    function initial_condition_transformed(x, t, aux_vars, equations)
        return Trixi.prim2cons(global2contravariant(initial_condition(x, t, equations),
                                                    aux_vars, equations), aux_vars,
                               equations)
    end
    return initial_condition_transformed
end

# Default version for standard Cartesian formulations
function transform_initial_condition(initial_condition, ::AbstractEquations)
    function initial_condition_transformed(x, t, equations)
        return Trixi.prim2cons(initial_condition(x, t, equations), equations)
    end
    return initial_condition_transformed
end

"""
    contravariant2global(u, aux_vars, equations)

Transform the vector `u` of solution variables with the momentum or velocity given in terms 
of local contravariant components into the global coordinate system specified by the 
`GlobalCoordinateSystem` type parameter in [`AbstractCovariantEquations`](@ref). `u` is a 
vector type of the correct length `nvariables(equations)`. Notice the function doesn't 
include any error checks for the purpose of efficiency, so please make sure your input is 
correct. The inverse conversion is performed by [`global2contravariant`](@ref).
"""
function contravariant2global end

"""
    global2contravariant(u, aux_vars, equations)

Transform the vector `u` of solution variables with momentum or velocity components
given with respect to the global coordinate system into local contravariant components. The 
global coordinate system is specified by the `GlobalCoordinateSystem` type parameter in 
[`AbstractCovariantEquations`](@ref). `u` is a vector type of the correct length 
`nvariables(equations)`. Notice the function doesn't include any error checks for the 
purpose of efficiency, so please make sure your input is correct. The inverse conversion is 
performed by [`contravariant2global`](@ref).
"""
function global2contravariant end

# By default, the equations are assumed to be formulated in Cartesian coordinates. This 
# function is specialized where needed.
function cartesian2global(u, x, equations::AbstractEquations)
    return u
end

# Default cons2cons and prim2cons methods which take in unused aux_vars variable
@inline Trixi.cons2cons(u, aux_vars, ::AbstractEquations) = u
@inline Trixi.prim2cons(u, aux_vars, ::AbstractEquations) = u

# For the covariant form, the auxiliary variables are the the NDIMS*NDIMS_AMBIENT entries 
# of the exact covariant basis matrix, the NDIMS*NDIMS_AMBIENT entries of the exact 
# contravariant basis matrix, the exact area element, the upper-triangular covariant and 
# contravariant metric tensor components, and the upper-triangular Christoffel symbols of 
# the second kind
@inline have_aux_node_vars(::AbstractCovariantEquations) = True()

# Add up the total number of auxiliary variables for equations in covariant form
@inline function n_aux_node_vars(::AbstractCovariantEquations{NDIMS,
                                                              NDIMS_AMBIENT}) where {
                                                                                     NDIMS,
                                                                                     NDIMS_AMBIENT
                                                                                     }
    nvars_basis_covariant = NDIMS_AMBIENT * NDIMS
    nvars_basis_contravariant = NDIMS * NDIMS_AMBIENT
    nvars_area_element = 1
    nvars_metric_covariant = NDIMS * (NDIMS + 1) ÷ 2
    nvars_metric_contravariant = NDIMS * (NDIMS + 1) ÷ 2
    nvars_christoffel = NDIMS * NDIMS * (NDIMS + 1) ÷ 2

    return nvars_basis_covariant +
           nvars_basis_contravariant +
           nvars_area_element +
           nvars_metric_covariant +
           nvars_metric_contravariant +
           nvars_christoffel
end

# Return auxiliary variable names for 2D covariant form
@inline function auxvarnames(::AbstractCovariantEquations{2})
    return ("basis_covariant[1,1]",         # e₁ ⋅ a₁
            "basis_covariant[2,1]",         # e₂ ⋅ a₁
            "basis_covariant[3,1]",         # e₃ ⋅ a₁
            "basis_covariant[1,2]",         # e₁ ⋅ a₂
            "basis_covariant[2,2]",         # e₂ ⋅ a₂
            "basis_covariant[3,2]",         # e₃ ⋅ a₂
            "basis_contravariant[1,1]",     # e₁ ⋅ a¹
            "basis_contravariant[2,1]",     # e₂ ⋅ a¹
            "basis_contravariant[3,1]",     # e₃ ⋅ a¹
            "basis_contravariant[1,2]",     # e₁ ⋅ a²
            "basis_contravariant[2,2]",     # e₂ ⋅ a²
            "basis_contravariant[3,2]",     # e₃ ⋅ a²
            "area_element",                 # J = √(G₁₁G₂₂ - G₁₂G₂₁) = ||a₁ × a₂||
            "metric_covariant[1,1]",        # G₁₁
            "metric_covariant[1,2]",        # G₁₂ = G₂₁
            "metric_covariant[2,2]",        # G₂₂
            "metric_contravariant[1,1]",    # G¹¹
            "metric_contravariant[1,2]",    # G¹² = G²¹
            "metric_contravariant[2,2]",    # G²²
            "christoffel_symbols[1][1,1]",  # Γ¹₁₁
            "christoffel_symbols[1][1,2]",  # Γ¹₁₂ = Γ¹₂₁
            "christoffel_symbols[1][2,2]",  # Γ¹₂₂
            "christoffel_symbols[2][1,1]",  # Γ²₁₁
            "christoffel_symbols[2][1,2]",  # Γ²₁₂ = Γ²₂₁
            "christoffel_symbols[2][2,2]")  # Γ²₂₂
end

# Extract the covariant basis vectors a_i from the auxiliary variables as a matrix A, 
# where A[:, i] contains the components of the ith covariant tangent basis vector with 
# respect to the global (Cartesian or spherical) coordinate system
@inline function basis_covariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{3, 2}(aux_vars[1], aux_vars[2], aux_vars[3],
                         aux_vars[4], aux_vars[5], aux_vars[6])
end

# Extract the contravariant basis vectors a^i from the auxiliary variables as a matrix B, 
# where B[i, :] contains the components of the ith contravariant tangent basis vector with 
# respect to the global (Cartesian or spherical) coordinate system
@inline function basis_contravariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{2, 3}(aux_vars[7], aux_vars[8],
                         aux_vars[9], aux_vars[10],
                         aux_vars[11], aux_vars[12])
end

# Extract the area element J = (det(AᵀA))^(1/2) from the auxiliary variables
@inline function area_element(aux_vars, ::AbstractCovariantEquations{2})
    return aux_vars[13]
end

# Extract the covariant metric tensor components Gᵢⱼ from the auxiliary variables
@inline function metric_covariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{2, 2}(aux_vars[14], aux_vars[15],
                         aux_vars[15], aux_vars[16])
end

# Extract the contravariant metric tensor components Gⁱʲ from the auxiliary variables
@inline function metric_contravariant(aux_vars, ::AbstractCovariantEquations{2})
    return SMatrix{2, 2}(aux_vars[17], aux_vars[18],
                         aux_vars[18], aux_vars[19])
end

# Extract the Christoffel symbols of the second kind Γⁱⱼₖ from the auxiliary variables
@inline function christoffel_symbols(aux_vars, ::AbstractCovariantEquations{2})
    return (SMatrix{2, 2}(aux_vars[20], aux_vars[21], aux_vars[21], aux_vars[22]),
            SMatrix{2, 2}(aux_vars[23], aux_vars[24], aux_vars[24], aux_vars[25]))
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

# Local Lax-Friedrichs dissipation for abstract covariant equations, where dissipation is 
# applied to all conservative variables and the wave speed may depend on auxiliary variables
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr, aux_vars_ll,
                                                              aux_vars_rr,
                                                              orientation_or_normal_direction,
                                                              equations::AbstractCovariantEquations)
    λ = dissipation.max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                  orientation_or_normal_direction, equations)
    return -0.5f0 * area_element(aux_vars_ll, equations) * λ * (u_rr - u_ll)
end

# Non-conservative two-point flux that returns a vector of zeros. This is used for systems 
# that have nonconservative terms when expressed in some formulations, but not others. For 
# example, to recover a standard weak formulation for CovariantShallowWaterEquations2D, one 
# must use volume_flux = (flux_central, flux_nonconservative_zeros). Once the bottom 
# topography source term for the covariant form is added, however, we will use that 
# instead, and this function will most likely no longer be needed.
@inline function flux_nonconservative_zeros(u_ll::SVector{NVARS, RealT},
                                            u_rr::SVector{NVARS, RealT},
                                            aux_vars_ll, aux_vars_rr,
                                            orientation_or_normal_direction,
                                            equations::AbstractCovariantEquations{2,
                                                                                    NDIMS_AMBIENT,
                                                                                    GlobalCoordinateSystem,
                                                                                    NVARS}) where {
                                                                                                     NDIMS_AMBIENT,
                                                                                                     GlobalCoordinateSystem,
                                                                                                     NVARS,
                                                                                                     RealT
                                                                                                     }
    return zeros(SVector{NVARS, RealT})
end

# Convert a vector from a global spherical to Cartesian basis representation. A tangent 
# vector will have vrad = 0.
@inline function spherical2cartesian(vlon, vlat, vrad, x)
    # compute longitude and latitude
    lon, lat = atan(x[2], x[1]), asin(x[3] / norm(x))

    # compute trigonometric functions
    sinlon, coslon = sincos(lon)
    sinlat, coslat = sincos(lat)

    # return Cartesian components
    vx = -sinlon * vlon - sinlat * coslon * vlat + coslat * coslon * vrad
    vy = coslon * vlon - sinlat * sinlon * vlat + coslat * sinlon * vrad
    vz = coslat * vlat + sinlat * vrad
    return vx, vy, vz
end

# Convert a vector from a global Cartesian to spherical basis representation. A tangent 
# vector will have vrad = 0.
@inline function cartesian2spherical(vx, vy, vz, x)
    # compute longitude and latitude
    lon, lat = atan(x[2], x[1]), asin(x[3] / norm(x))

    # compute trigonometric functions
    sinlon, coslon = sincos(lon)
    sinlat, coslat = sincos(lat)

    # return spherical components
    vlon = -sinlon * vx + coslon * vy
    vlat = -sinlat * coslon * vx - sinlat * sinlon * vy + coslat * vz
    vrad = coslat * coslon * vx + coslat * sinlon * vy + sinlat * vz

    return vlon, vlat, vrad
end

abstract type AbstractCompressibleMoistEulerEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end

include("covariant_advection.jl")
include("covariant_shallow_water.jl")
include("compressible_moist_euler_2d_lucas.jl")
include("shallow_water_3d.jl")
include("reference_data.jl")
end # @muladd
