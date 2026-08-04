@muladd begin
#! format: noindent

@doc raw"""
    CovariantEulerEquations3D{GlobalCoordinateSystem} <:  
        AbstractCovariantEquations{3, 3, GlobalCoordinateSystem, 5}

Denoting the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) by 
$\nabla_j$ and summing over repeated indices, the compressible Euler equations with gravity 
can be expressed in covariant form as
```math
\begin{aligned}
\partial_t \rho + \nabla_j (\rho v^j) &= 0,\\
\partial_t (\rho v^i) + \nabla_j \tau^{ij} + \rho G^{ij}\partial_j \Phi &= 0,\\
\partial_t E + \nabla_j\big(v^j(E + p)\big) &= 0,
\end{aligned}
```
where $\rho$ is the density, $v^i$ and $G^{ij}$ are the contravariant velocity and metric
tensor components, $\Phi$ is the geopotential, and $\partial_j$ is used as a shorthand for
$\partial / \partial \xi^j$. The total energy density $E$ is taken to include the
gravitational potential energy, such that the pressure $p$ is given in terms of the
conservative variables by the equation of state
```math
p = (\gamma - 1)\left(E - \frac{1}{2}\rho G_{ij} v^i v^j - \rho \Phi\right),
```
with $\gamma$ denoting the ratio of specific heats. Because the potential energy is already 
accounted for within $E$, no additional source term due to gravity appears in the energy 
equation. The contravariant momentum flux tensor components are given by
```math
\tau^{ij} = \rho v^i v^j + p G^{ij}.
```
As with the covariant shallow water equations (see 
[`CovariantShallowWaterEquations2D`](@ref)), this system may be formulated on the reference 
element as a system of conservation laws with a source term, as given by
```math
J \frac{\partial}{\partial t}
\left[\begin{array}{c} \rho \\ \rho v^1 \\ \rho v^2 \\ E \end{array}\right] 
+
\frac{\partial}{\partial \xi^1} 
\left[\begin{array}{c} J \rho v^1 \\ J \tau^{11} \\ J \tau^{12} \\ J v^1 (E + p) \end{array}\right]
+ 
\frac{\partial}{\partial \xi^2} 
\left[\begin{array}{c} J \rho v^2 \\ J \tau^{21} \\ J \tau^{22} \\ J v^2 (E + p) \end{array}\right] 
= J \left[\begin{array}{c} 0 \\ 
-\Gamma^1_{jk}\tau^{jk} - J\rho G^{1j}\partial_j \Phi \\ 
-\Gamma^2_{jk}\tau^{jk} - J\rho G^{2j}\partial_j \Phi \\
0
 \end{array}\right],
```
where the Christoffel symbols of the second kind $\Gamma^i_{jk}$ are defined as in 
[`CovariantShallowWaterEquations2D`](@ref).

!!! note
    The geometric part of the source term above, involving the Christoffel symbols
    $\Gamma^i_{jk}$, is **not currently implemented**. Only the gravitational contribution 
    to the momentum source term is included, through the exported function 
    `source_terms_gravity`.

## References
- M. Baldauf and F. Prill. Comparison of two Euler equation sets in a Discontinuous
  Galerkin solver for atmospheric modelling (BRIDGE v0.9).
"""
struct CovariantEulerEnergyEquations3D{GlobalCoordinateSystem,
                                          RealT <: Real} <:
       AbstractCovariantEulerEquations{3, GlobalCoordinateSystem, 5}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    gravity::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
    global_coordinate_system::GlobalCoordinateSystem
    function CovariantEulerEnergyEquations3D(; c_p, c_v,
                                             gravity,
                                             global_coordinate_system = GlobalCartesianCoordinates())
        c_p, c_v, gravity = promote(c_p, c_v, gravity)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        K = p_0 * (R / p_0)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        return new{typeof(global_coordinate_system), typeof(c_p)}(p_0, c_p, c_v, gravity, R,
                                                                  gamma,
                                                                  inv_gamma_minus_one,
                                                                  K, stolarsky_factor)
    end
end

have_nonconservative_terms(::CovariantEulerEnergyEquations3D) = False()

# The conservative variables are the height and contravariant momentum components
function varnames(::typeof(cons2cons), ::CovariantEulerEnergyEquations3D)
    return ("rho", "rho_vcon1", "rho_vcon2", "rho_vcon3", "rho_e")
end

# The primitive variables are the height and contravariant velocity components
function varnames(::typeof(cons2prim), ::CovariantEulerEnergyEquations3D)
    return ("rho", "vcon1", "vcon2", "vcon3", "p")
end

function varnames(::typeof(contravariant2global),
                  ::CovariantEulerEnergyEquations3D)
    return ("rho", "v1", "v2", "v3", "rho_e")
end

# Convenience functions to extract physical variables from state vector
@inline density(u, ::CovariantEulerEnergyEquations3D) = u[1]

@inline velocity_contravariant(u, ::CovariantEulerEnergyEquations3D) = SVector(u[2] /
                                                                               u[1],
                                                                               u[3] /
                                                                               u[1],
                                                                               u[4] /
                                                                               u[1])

@inline momentum_contravariant(u, ::CovariantEulerEnergyEquations3D) = SVector(u[2],
                                                                               u[3],
                                                                               u[4])

@inline total_energy(u, ::CovariantEulerEnergyEquations3D) = u[5]

@inline energy_density(u, ::CovariantEulerEnergyEquations3D) = u[5] / u[1]

@inline function cons2prim(u, aux_vars,
                           equations::CovariantEulerEnergyEquations3D)
    rho = density(u, equations)
    vcon = velocity_contravariant(u, equations)
    p = pressure(u, aux_vars, equations)
    return SVector(rho, vcon[1], vcon[2], vcon[3], p)
end

@inline function prim2cons(u, aux_vars,
                           equations::CovariantEulerEnergyEquations3D)
    rho, vcon1, vcon2, vcon3, p = u
    vcon = SVector(vcon1, vcon2, vcon3)
    Gcov = metric_covariant(aux_vars, equations)
    ekin = 0.5f0 * dot(Gcov * vcon, vcon) * rho
    phi = geopotential(aux_vars, equations)
    rho_e_total = p * equations.inv_gamma_minus_one + ekin + rho * phi
    return SVector(rho, rho * vcon1, rho * vcon2, rho * vcon3, rho_e_total)
end

@inline function cons2entropy(u, aux_vars,
                              equations::CovariantEulerEnergyEquations3D)
    Gcov = metric_covariant(aux_vars, equations)
    rho, rho_vcon1, rho_vcon2, rho_vcon3, rho_e_total = u
    p = pressure(u, aux_vars, equations)
    vcon = SVector(rho_vcon1 / rho, rho_vcon2 / rho, rho_vcon3 / rho)
    vcov = Gcov * vcon
    s = log(p / rho^equations.gamma)

    # Entropy variables (covariant form)
    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         0.5f0 * rho / p * dot(vcov, vcon)
    w2 = rho / p * vcov[1]
    w3 = rho / p * vcov[2]
    w4 = rho / p * vcov[3]
    w5 = -rho / p

    return SVector(w1, w2, w3, w4, w5)
end

# Convert contravariant momentum components to the global coordinate system
@inline function contravariant2global(u, aux_vars,
                                      equations::CovariantEulerEnergyEquations3D)
    v1, v2, v3 = basis_covariant(aux_vars, equations) *
                 SVector(u[2], u[3], u[4])
    return SVector(u[1], v1, v2, v3, u[5])
end

# Convert momentum components in the global coordinate system to contravariant components
@inline function global2contravariant(u, aux_vars,
                                      equations::CovariantEulerEnergyEquations3D)
    vcon1, vcon2, vcon3 = basis_contravariant(aux_vars, equations) *
                   SVector(u[2], u[3], u[4])
    return SVector(u[1], vcon1, vcon2, vcon3, u[5])
end

# Maximum wave speeds with respect to the covariant basis
@inline function max_abs_speeds(u, aux_vars,
                                equations::CovariantEulerEnergyEquations3D)
    rho, vcon1, vcon2, vcon3, p = cons2prim(u, aux_vars, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    c1 = sqrt(equations.gamma * p / rho * Gcon[1, 1])
    c2 = sqrt(equations.gamma * p / rho * Gcon[2, 2])
    c3 = sqrt(equations.gamma * p / rho * Gcon[3, 3])

    return abs(vcon1) + c1, abs(vcon2) + c2, abs(vcon3) + c3
end

# Maximum wave speed along the normal direction in reference space
@inline function max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               normal_direction::AbstractVector,
                               equations::CovariantEulerEnergyEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_vars_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_vars_rr, equations)

    # Calculate the normal velocities and sound speeds
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2]
            +
            v3_ll * normal_direction[3])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll *
                (normal_direction' * metric_contravariant(aux_vars_ll, equations) *
                 normal_direction))

    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2]
            +
            v3_rr * normal_direction[3])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr *
                (normal_direction' * metric_contravariant(aux_vars_rr, equations) *
                 normal_direction))

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

@inline function boundary_condition_slip_wall(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function,
                                              equations::CovariantEulerEnergyEquations3D)
    # Reflect the normal contravariant velocity component and keep the tangential component unchanged
    u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], -u_inner[4], u_inner[5])

    # Compute the flux at the boundary using the surface flux function
    flux = surface_flux_function(u_inner, u_boundary, aux_inner, aux_inner,
                                 normal_direction, equations)
    return flux
end

# Standard geometric and Coriolis source terms for a rotating sphere
@inline function source_terms_geometric_coriolis(u, x, t, aux_vars,
                                                 equations::CovariantEulerEnergyEquations3D)
    # Geometric variables
    Gcon = metric_contravariant(aux_vars, equations)
    Gamma1, Gamma2, Gamma3 = christoffel_symbols(aux_vars, equations)
    J = volume_element(aux_vars, equations)

    # Physical variables
    rho = density(u, equations)
    Omega = basis_contravariant(aux_vars, equations)[:, 3]
    rho_vcon = momentum_contravariant(u, equations)
    vcon = velocity_contravariant(u, equations)
    p = pressure(u, aux_vars, equations)

    # Doubly-contravariant momentum flux tensor
    momentum_flux = rho_vcon * vcon' + p * Gcon

    # Geometric source term
    s_geo = SVector(sum(Gamma1 .* momentum_flux), sum(Gamma2 .* momentum_flux), sum(Gamma3 .* momentum_flux))

    # Coriolis source term
    cross_product = SVector(Omega[2] * rho_vcon[3] - Omega[3] * rho_vcon[2],
                            Omega[3] * rho_vcon[1] - Omega[1] * rho_vcon[3],
                            Omega[1] * rho_vcon[2] - Omega[2] * rho_vcon[1])
    s_cor = 2 * EARTH_ROTATION_RATE * J * Gcon * cross_product

    # Bouyancy source term
    r = norm(x)
    GM = EARTH_RADIUS^2 * EARTH_GRAVITATIONAL_ACCELERATION
    Gcov = metric_covariant(aux_vars, equations)
    derivative_phi = -GM / r^2
    s_boy = [0, 0, derivative_phi * rho * Gcon[3, 3] * sqrt(Gcov[3, 3])]

    # Combined source terms
    source_1 = s_geo[1] + s_cor[1]
    source_2 = s_geo[2] + s_cor[2]
    source_3 = s_geo[3] + s_cor[3] + s_boy[3]

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(0, -source_1, -source_2, -source_3, 0)
end
end # @muladd
