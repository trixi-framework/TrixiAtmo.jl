@muladd begin
#! format: noindent

@doc raw"""
    CovariantShallowWaterEquations2D{GlobalCoordinateSystem} <:  
        AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}

Denoting the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) by 
$\nabla_j$ and summing over repeated indices, the shallow water equations can be expressed 
on a two-dimensional surface in three-dimensional ambient space as
```math
\begin{aligned}
\partial_t h + \nabla_j (hv^j) &= 0,\\
\partial_t (hv^i) + \nabla_j \tau^{ij} + gh G^{ij}\partial_j b
&= -fJ G^{ij}\varepsilon_{jk} hv^k,
\end{aligned}
```
where $h$ is the geopotential height (equal to the total geopotential height $H$ for zero
bottom topography), $v^i$ and $G^{ij}$ are the contravariant velocity and metric tensor
components, $g$ is the gravitational acceleration, $b$ is the bottom topography, $f$ is the 
Coriolis parameter, $J$ is the area element, $\varepsilon$ is the Levi-Civita symbol, and 
$\partial_j$ is used as a shorthand for $\partial / \partial \xi^j$. The contravariant 
momentum flux tensor components are given by
```math
\tau^{ij} = hv^i v^j + \frac{1}{2}G^{ij}gh^2.
```
The covariant shallow water equations with constant bottom topography can be formulated on 
the reference element as a system of conservation laws with a source term (implemented in 
the exported function `source_terms_geometric_coriolis`), as given by
```math
J \frac{\partial}{\partial t}
\left[\begin{array}{c} h \\ hv^1 \\ hv^2 \end{array}\right] 
+
\frac{\partial}{\partial \xi^1} 
\left[\begin{array}{c} J h v^1 \\ J \tau^{11} \\ J \tau^{12} \end{array}\right]
+ 
\frac{\partial}{\partial \xi^2} 
\left[\begin{array}{c} J h v^2 \\ J \tau^{21} \\ J \tau^{22}  \end{array}\right] 
= J \left[\begin{array}{c} 0 \\ 
-\Gamma^1_{jk}\tau^{jk} - f J \big(G^{12}hv^1 - G^{11}hv^2\big) \\ 
-\Gamma^2_{jk}\tau^{jk} - f J \big(G^{22}hv^1 - G^{21}hv^2\big)
 \end{array}\right].
```
Note that the geometric contribution to the source term involves the Christoffel symbols of
the second kind, which can been expressed in terms of the covariant metric tensor 
components $G_{ij}$ as 
```math
\Gamma_{jk}^i = 
\frac{1}{2}G^{il}\big(\partial_j G_{kl} + \partial_k G_{jl} - \partial_l G_{jk}\big).
```
## References
- M. Baldauf (2020). Discontinuous Galerkin solver for the shallow-water equations in
  covariant form on the sphere and the ellipsoid. Journal of Computational Physics 
  410:109384. [DOI: 10.1016/j.jcp.2020.109384](https://doi.org/10.1016/j.jcp.2020.109384) 
- L. Bao, R. D. Nair, and H. M. Tufo (2014). A mass and momentum flux-form high-order
  discontinuous Galerkin shallow water model on the cubed-sphere. A mass and momentum 
  flux-form high-order discontinuous Galerkin shallow water model on the cubed-sphere. 
  Journal of Computational Physics 271:224-243. 
  [DOI: 10.1016/j.jcp.2013.11.033](https://doi.org/10.1016/j.jcp.2013.11.033)

!!! note
    When solving problems with variable bottom topography as well as when using
    entropy-stable schemes, [SplitCovariantShallowWaterEquations2D](@ref) should be used
    instead.
"""
struct CovariantShallowWaterEquations2D{GlobalCoordinateSystem, RealT <: Real} <:
       AbstractCovariantShallowWaterEquations2D{GlobalCoordinateSystem}
    gravity::RealT  # acceleration due to gravity
    rotation_rate::RealT  # rotation rate for Coriolis term 
    global_coordinate_system::GlobalCoordinateSystem
    function CovariantShallowWaterEquations2D(gravity::RealT,
                                              rotation_rate::RealT;
                                              global_coordinate_system = GlobalCartesianCoordinates()) where {RealT <:
                                                                                                              Real}
        return new{typeof(global_coordinate_system), RealT}(gravity, rotation_rate)
    end
end

# Until we implement bottom topography, there are no nonconservative terms
Trixi.have_nonconservative_terms(::CovariantShallowWaterEquations2D) = False()

# The conservative variables are the height and contravariant momentum components
function Trixi.varnames(::typeof(cons2cons), ::AbstractCovariantShallowWaterEquations2D)
    return ("h", "h_vcon1", "h_vcon2")
end

# The primitive variables are the height and contravariant velocity components
function Trixi.varnames(::typeof(cons2prim), ::AbstractCovariantShallowWaterEquations2D)
    return ("H", "vcon1", "vcon2")
end

# The change of variables contravariant2global converts the two local contravariant vector 
# components u[2] and u[3] to the three global vector components specified by 
# equations.global_coordinate_system (e.g. spherical or Cartesian). This transformation 
# works for both primitive and conservative variables, although varnames refers 
# specifically to transformations from conservative variables.
function Trixi.varnames(::typeof(contravariant2global),
                        ::AbstractCovariantShallowWaterEquations2D)
    return ("h", "h_v1", "h_v2", "h_v3")
end

# Convenience functions to extract physical variables from state vector
@inline Trixi.waterheight(u, ::AbstractCovariantShallowWaterEquations2D) = u[1]
@inline velocity_contravariant(u,
::AbstractCovariantShallowWaterEquations2D) = SVector(u[2] /
                                                      u[1],
                                                      u[3] /
                                                      u[1])
@inline momentum_contravariant(u,
::AbstractCovariantShallowWaterEquations2D) = SVector(u[2],
                                                      u[3])

@inline function Trixi.cons2prim(u, aux_vars,
                                 equations::AbstractCovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u
    h_s = bottom_topography(aux_vars, equations)
    return SVector(h + h_s, h_vcon1 / h, h_vcon2 / h)
end

@inline function Trixi.prim2cons(u, aux_vars,
                                 equations::AbstractCovariantShallowWaterEquations2D)
    H, vcon1, vcon2 = u
    h_s = bottom_topography(aux_vars, equations)
    h = H - h_s
    return SVector(h, h * vcon1, h * vcon2)
end

# Entropy variables are w = (g(h+hₛ) - (v₁v¹ + v₂v²)/2, v₁, v₂)ᵀ
@inline function Trixi.cons2entropy(u, aux_vars,
                                    equations::AbstractCovariantShallowWaterEquations2D)
    h = waterheight(u, equations)
    h_s = bottom_topography(aux_vars, equations)
    vcon = velocity_contravariant(u, equations)
    vcov = metric_covariant(aux_vars, equations) * vcon
    return SVector{3}(equations.gravity * (h + h_s) - 0.5f0 * dot(vcov, vcon),
                      vcov[1], vcov[2])
end

# Convert contravariant momentum components to the global coordinate system
@inline function contravariant2global(u, aux_vars,
                                      equations::AbstractCovariantShallowWaterEquations2D)
    h_v1, h_v2, h_v3 = basis_covariant(aux_vars, equations) *
                       momentum_contravariant(u, equations)
    return SVector(waterheight(u, equations), h_v1, h_v2, h_v3)
end

# Convert momentum components in the global coordinate system to contravariant components
@inline function global2contravariant(u, aux_vars,
                                      equations::AbstractCovariantShallowWaterEquations2D)
    h_vcon1, h_vcon2 = basis_contravariant(aux_vars, equations) *
                       SVector(u[2], u[3], u[4])
    return SVector(u[1], h_vcon1, h_vcon2)
end

# Entropy function (total energy) given by S = (h(v₁v¹ + v₂v²) + gh² + ghhₛ)/2
@inline function Trixi.entropy(u, aux_vars,
                               equations::AbstractCovariantShallowWaterEquations2D)
    h = waterheight(u, equations)
    h_s = bottom_topography(aux_vars, equations)
    vcon = velocity_contravariant(u, equations)
    vcov = metric_covariant(aux_vars, equations) * vcon
    return 0.5f0 * (h * dot(vcov, vcon) + equations.gravity * h^2) +
           equations.gravity * h * h_s
end

# Flux as a function of the state vector u, as well as the auxiliary variables aux_vars, 
# which contain the geometric information required for the covariant form
@inline function Trixi.flux(u, aux_vars, orientation::Integer,
                            equations::CovariantShallowWaterEquations2D)
    # Geometric variables
    Gcon = metric_contravariant(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    h = waterheight(u, equations)
    h_vcon = momentum_contravariant(u, equations)

    # Compute and store the velocity and gravitational terms
    vcon = h_vcon[orientation] / h
    gravitational_term = 0.5f0 * equations.gravity * h^2

    # Compute the momentum flux components in the desired orientation
    momentum_flux_1 = h_vcon[1] * vcon + Gcon[1, orientation] * gravitational_term
    momentum_flux_2 = h_vcon[2] * vcon + Gcon[2, orientation] * gravitational_term

    return SVector(J * h_vcon[orientation], J * momentum_flux_1, J * momentum_flux_2)
end

# Standard geometric and Coriolis source terms for a rotating sphere
@inline function source_terms_geometric_coriolis(u, x, t, aux_vars,
                                                 equations::CovariantShallowWaterEquations2D)
    # Geometric variables
    Gcon = metric_contravariant(aux_vars, equations)
    Gamma1, Gamma2 = christoffel_symbols(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    h = waterheight(u, equations)
    h_vcon = momentum_contravariant(u, equations)
    v_con = velocity_contravariant(u, equations)

    # Doubly-contravariant flux tensor
    momentum_flux = h_vcon * v_con' + 0.5f0 * equations.gravity * h^2 * Gcon

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Geometric source term
    s_geo = SVector(sum(Gamma1 .* momentum_flux), sum(Gamma2 .* momentum_flux))

    # Combined source terms
    source_1 = s_geo[1] + f * J * (Gcon[1, 2] * h_vcon[1] - Gcon[1, 1] * h_vcon[2])
    source_2 = s_geo[2] + f * J * (Gcon[2, 2] * h_vcon[1] - Gcon[2, 1] * h_vcon[2])

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(zero(eltype(u)), -source_1, -source_2)
end

# Maximum wave speed along the normal direction in reference space
@inline function Trixi.max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                     orientation,
                                     equations::AbstractCovariantShallowWaterEquations2D)
    # Geometric variables
    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    Gcon_rr = metric_contravariant(aux_vars_rr, equations)

    # Physical variables 
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    h_vcon_ll = momentum_contravariant(u_ll, equations)
    h_vcon_rr = momentum_contravariant(u_rr, equations)

    return max(abs(h_vcon_ll[orientation] / h_ll) +
               sqrt(Gcon_ll[orientation, orientation] * h_ll * equations.gravity),
               abs(h_vcon_rr[orientation] / h_rr) +
               sqrt(Gcon_rr[orientation, orientation] * h_rr * equations.gravity))
end

# Maximum wave speeds with respect to the covariant basis
@inline function Trixi.max_abs_speeds(u, aux_vars,
                                      equations::AbstractCovariantShallowWaterEquations2D)
    vcon = velocity_contravariant(u, equations)
    h = waterheight(u, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    return abs(vcon[1]) + sqrt(Gcon[1, 1] * h * equations.gravity),
           abs(vcon[2]) + sqrt(Gcon[2, 2] * h * equations.gravity)
end

# If the initial velocity field is defined in Cartesian coordinates and the chosen global 
# coordinate system is spherical, perform the appropriate conversion
@inline function cartesian2global(u, x,
                                  ::AbstractCovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
    h_vlon, h_vlat, h_vrad = cartesian2spherical(u[2], u[3], u[4], x)
    return SVector(u[1], h_vlon, h_vlat, h_vrad)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is Cartesian, perform the appropriate conversion
@inline function spherical2global(u, x,
                                  ::AbstractCovariantShallowWaterEquations2D{GlobalCartesianCoordinates})
    h_vx, h_vy, h_vz = spherical2cartesian(u[2], u[3], u[4], x)
    return SVector(u[1], h_vx, h_vy, h_vz)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is spherical, do not convert
@inline function spherical2global(u, x,
                                  ::AbstractCovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
    return u
end
end # @muladd
