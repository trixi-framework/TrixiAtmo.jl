@muladd begin
#! format: noindent

@doc raw"""
    CovariantShallowWaterEquations2D{GlobalCoordinateSystem} <:  
        AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}

Denoting the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) by 
$\nabla_b$ and summing over repeated indices, the shallow water equations can be expressed 
on a two-dimensional surface in three-dimensional ambient space as
```math
\begin{aligned}
\partial_t h + \nabla_b (hv^b) &= 0,\\
\partial_t (hv^a) + \nabla_b (hv^av^b) + gh G^{ab}\partial_b(h + b) 
&= -fJ G^{ab}\varepsilon_{bc} hv^c,
\end{aligned}
```
where $h$ is the fluid height, $v^a$ and $G^{ab}$ are the contravariant velocity and metric 
tensor components, $g$ is the gravitational constant, $f$ is the Coriolis parameter, and 
$J$ is the area element. Combining the inertial and gravitational terms in order to define 
the momentum flux components
```math
\tau^{ab} = hv^a v^b + \frac{1}{2}G^{ab}gh^2,
```
the covariant shallow water equations can be expressed as a system of conservation laws 
with a source term:
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
-\Gamma^1_{ab}\tau^{ab} - f J \big(G^{12}hv^1 - G^{11}hv^2\big) \\ 
-\Gamma^2_{ab}\tau^{ab} - f J \big(G^{22}hv^1 - G^{21}hv^2\big)
 \end{array}\right].
```
TrixiAtmo.jl implements standard weak-form DG methods for both the above system, as well as 
entropy-stable split forms based on a novel flux-differencing discretization of the 
covariant derivative.

## References
- M. Baldauf (2020). Discontinuous Galerkin solver for the shallow-water equations in
  covariant form on the sphere and the ellipsoid. Journal of Computational Physics 
  410:109384. [DOI: 10.1016/j.jcp.2020.109384](https://doi.org/10.1016/j.jcp.2020.109384) 
- L. Bao, R. D. Nair, and H. M. Tufo (2014). A mass and momentum flux-form high-order
  discontinuous Galerkin shallow water model on the cubed-sphere. A mass and momentum 
  flux-form high-order discontinuous Galerkin shallow water model on the cubed-sphere. 
  Journal of Computational Physics 271:224-243. 
  [DOI: 10.1016/j.jcp.2013.11.033](https://doi.org/10.1016/j.jcp.2013.11.033)

!!! warning "Experimental implementation"
    The use of entropy-stable split-form/flux-differencing formulations for covariant 
    equations is an experimental feature and may change in future releases.
"""
struct CovariantShallowWaterEquations2D{GlobalCoordinateSystem, RealT <: Real} <:
       AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}
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

# Our implementation of flux-differencing formulation uses nonconservative terms, but the 
# standard weak form does not. To handle both options, we have defined a dummy kernel for 
# the nonconservative terms that does nothing when VolumeIntegralWeakForm is used with a 
# nonconservative system.
Trixi.have_nonconservative_terms(::CovariantShallowWaterEquations2D) = True()

# The conservative variables are the height and contravariant momentum components
function Trixi.varnames(::typeof(cons2cons), ::CovariantShallowWaterEquations2D)
    return ("h", "h_vcon1", "h_vcon2")
end

# The primitive variables are the height and contravariant velocity components
function Trixi.varnames(::typeof(cons2prim), ::CovariantShallowWaterEquations2D)
    return ("h", "vcon1", "vcon2")
end

# The change of variables contravariant2global converts the two local contravariant vector 
# components u[2] and u[3] to the three global vector components specified by 
# equations.global_coordinate_system (e.g. spherical or Cartesian). This transformation 
# works for both primitive and conservative variables, although varnames refers 
# specifically to transformations from conservative variables.
function Trixi.varnames(::typeof(contravariant2global),
                        ::CovariantShallowWaterEquations2D)
    return ("h", "h_vglo1", "h_vglo2", "h_vglo3")
end

# Convenience functions to extract physical variables from state vector
@inline waterheight(u, ::CovariantShallowWaterEquations2D) = u[1]
@inline velocity_contravariant(u, ::CovariantShallowWaterEquations2D) = SVector(u[2] /
                                                                                u[1],
                                                                                u[3] /
                                                                                u[1])
@inline momentum_contravariant(u, ::CovariantShallowWaterEquations2D) = SVector(u[2],
                                                                                u[3])

@inline function Trixi.cons2prim(u, aux_vars, ::CovariantShallowWaterEquations2D)
    h, h_vcon1, h_vcon2 = u
    return SVector(h, h_vcon1 / h, h_vcon2 / h)
end

@inline function Trixi.prim2cons(u, aux_vars, ::CovariantShallowWaterEquations2D)
    h, vcon1, vcon2 = u
    return SVector(h, h * vcon1, h * vcon2)
end

@inline function Trixi.cons2entropy(u, aux_vars,
                                    equations::CovariantShallowWaterEquations2D)
    h = waterheight(u, equations)
    vcon = velocity_contravariant(u, equations)
    vcov = metric_covariant(aux_vars, equations) * vcon
    return SVector{3}(equations.gravity * h - 0.5f0 * dot(vcov, vcon), vcov[1], vcov[2])
end

# Convert contravariant momentum components to the global coordinate system
@inline function contravariant2global(u, aux_vars,
                                      equations::CovariantShallowWaterEquations2D)
    h_vglo1, h_vglo2, h_vglo3 = basis_covariant(aux_vars, equations) *
                                momentum_contravariant(u, equations)
    return SVector(waterheight(u, equations), h_vglo1, h_vglo2, h_vglo3)
end

# Convert momentum components in the global coordinate system to contravariant components
@inline function global2contravariant(u, aux_vars,
                                      equations::CovariantShallowWaterEquations2D)
    h_vcon1, h_vcon2 = basis_contravariant(aux_vars, equations) *
                       SVector(u[2], u[3], u[4])
    return SVector(u[1], h_vcon1, h_vcon2)
end

# Entropy function (total energy per unit volume)
@inline function Trixi.entropy(u, aux_vars, equations::CovariantShallowWaterEquations2D)
    h = waterheight(u, equations)
    vcon = velocity_contravariant(u, equations)
    vcov = metric_covariant(aux_vars, equations) * vcon
    return 0.5f0 * (dot(vcov, vcon) + equations.gravity * h^2)
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

# Symmetric part of entropy-conservative flux
@inline function Trixi.flux_ec(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               orientation::Integer,
                               equations::CovariantShallowWaterEquations2D)
    # Geometric variables
    J_ll = area_element(aux_vars_ll, equations)
    J_rr = area_element(aux_vars_rr, equations)

    # Physical variables
    h_vcon_ll = momentum_contravariant(u_ll, equations)
    h_vcon_rr = momentum_contravariant(u_rr, equations)
    vcon_ll = velocity_contravariant(u_ll, equations)
    vcon_rr = velocity_contravariant(u_rr, equations)

    # Scaled mass flux in conservative form
    mass_flux = 0.5f0 * (J_ll * h_vcon_ll[orientation] + J_rr * h_vcon_rr[orientation])

    # Half of scaled inertial flux in conservative form
    momentum_flux = 0.25f0 * (J_ll * h_vcon_ll * vcon_ll[orientation] +
                     J_rr * h_vcon_rr * vcon_rr[orientation])

    return SVector(mass_flux, momentum_flux[1], momentum_flux[2])
end

# Entropy-conservative flux with local Lax-Friedrichs dissipation
const flux_ec_llf = FluxPlusDissipation(flux_ec,
                                        DissipationLocalLaxFriedrichs(max_abs_speed_naive))

# Non-symmetric part of entropy-conservative flux
@inline function flux_nonconservative_ec(u_ll, u_rr, aux_vars_ll,
                                         aux_vars_rr,
                                         orientation::Integer,
                                         equations::CovariantShallowWaterEquations2D)
    # Geometric variables
    Gcov_ll = metric_covariant(aux_vars_ll, equations)
    Gcov_rr = metric_covariant(aux_vars_rr, equations)
    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    J_ll = area_element(aux_vars_ll, equations)
    J_rr = area_element(aux_vars_rr, equations)

    # Physical variables
    h_vcon_ll = momentum_contravariant(u_ll, equations)
    h_vcon_rr = momentum_contravariant(u_rr, equations)
    vcov_ll = Gcov_ll * velocity_contravariant(u_ll, equations)
    vcov_rr = Gcov_rr * velocity_contravariant(u_rr, equations)

    # Half of inertial term in non-conservative form
    nonconservative_term_inertial = 0.5f0 * Gcon_ll *
                                    (J_ll * h_vcon_ll[orientation] * vcov_rr +
                                     J_rr * h_vcon_rr[orientation] * vcov_ll)

    # Gravity term in non-conservative form
    nonconservative_term_gravitational = equations.gravity * J_ll *
                                         Gcon_ll[:, orientation] *
                                         waterheight(u_ll, equations) *
                                         waterheight(u_rr, equations)

    nonconservative_term = nonconservative_term_inertial +
                           nonconservative_term_gravitational

    return SVector(zero(eltype(u_ll)), nonconservative_term[1], nonconservative_term[2])
end

# Geometric and Coriolis sources for a rotating sphere with VolumeIntegralWeakForm
@inline function source_terms_weak_form(u, x, t, aux_vars,
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

# Geometric and Coriolis sources for a rotating sphere with VolumeIntegralFluxDifferencing
@inline function source_terms_ec(u, x, t, aux_vars,
                                 equations::CovariantShallowWaterEquations2D)
    # Geometric variables
    Gcov = metric_covariant(aux_vars, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    (Gamma1, Gamma2) = christoffel_symbols(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    h_vcon = momentum_contravariant(u, equations)
    v_con = velocity_contravariant(u, equations)

    # Doubly-contravariant and mixed inertial flux tensors
    h_vcon_vcon = h_vcon * v_con'
    h_vcov_vcon = Gcov * h_vcon_vcon

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Geometric source term
    s_geo = 0.5f0 * (SVector(sum(Gamma1 .* h_vcon_vcon), sum(Gamma2 .* h_vcon_vcon)) -
             Gcon * (Gamma1 * h_vcov_vcon[1, :] + Gamma2 * h_vcov_vcon[2, :]))

    # Combined source terms
    source_1 = s_geo[1] + f * J * (Gcon[1, 2] * h_vcon[1] - Gcon[1, 1] * h_vcon[2])
    source_2 = s_geo[2] + f * J * (Gcon[2, 2] * h_vcon[1] - Gcon[2, 1] * h_vcon[2])

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(zero(eltype(u)), -source_1, -source_2)
end

# Maximum wave speed along the normal direction in reference space
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                           orientation,
                                           equations::CovariantShallowWaterEquations2D)
    # Geometric variables
    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    Gcon_rr = metric_contravariant(aux_vars_rr, equations)

    # Physical variables 
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_ll, equations)
    h_vcon_ll = momentum_contravariant(u_ll, equations)
    h_vcon_rr = momentum_contravariant(u_rr, equations)

    return max(abs(h_vcon_ll[orientation] / h_ll) +
               sqrt(Gcon_ll[orientation, orientation] * h_ll * equations.gravity),
               abs(h_vcon_rr[orientation] / h_rr) +
               sqrt(Gcon_rr[orientation, orientation] * h_rr * equations.gravity))
end

# Maximum wave speeds with respect to the covariant basis
@inline function Trixi.max_abs_speeds(u, aux_vars,
                                      equations::CovariantShallowWaterEquations2D)
    vcon = velocity_contravariant(u, equations)
    h = waterheight(u, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    return abs(vcon[1]) + sqrt(Gcon[1, 1] * h * equations.gravity),
           abs(vcon[2]) + sqrt(Gcon[2, 2] * h * equations.gravity)
end

# If the initial velocity field is defined in Cartesian coordinates and the chosen global 
# coordinate system is spherical, perform the appropriate conversion
@inline function cartesian2global(u, x,
                                  ::CovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
    h_vlon, h_vlat, h_vrad = cartesian2spherical(u[2], u[3], u[4], x)
    return SVector(u[1], h_vlon, h_vlat, h_vrad)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is Cartesian, perform the appropriate conversion
@inline function spherical2global(u, x,
                                  ::CovariantShallowWaterEquations2D{GlobalCartesianCoordinates})
    h_vx, h_vy, h_vz = spherical2cartesian(u[2], u[3], u[4], x)
    return SVector(u[1], h_vx, h_vy, h_vz)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is spherical, do not convert
@inline function spherical2global(u, x,
                                  ::CovariantShallowWaterEquations2D{GlobalSphericalCoordinates})
    return u
end
end # @muladd
