@muladd begin
#! format: noindent

@doc raw"""
    CompressibleEulerVectorInvariantEquationsWithGravity3D(; c_p, c_v, gravity)

The compressible Euler equations with gravity in vector-invariant form with the potential
temperature as prognostic variable
```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{v}) &= 0, \\
\frac{\partial \vec{v}}{\partial t} + \nabla \left( \frac{1}{2} \lVert \vec{v} \rVert^2 \right)
+ (\nabla \times \vec{v}) \times \vec{v} + c_p \theta \nabla \pi + \nabla \phi &= 0, \\
\frac{\partial (\rho \theta)}{\partial t} + \nabla \cdot (\rho \theta \vec{v}) &= 0,
\end{aligned}
```
for an ideal gas with ratio of specific heats ``\gamma = c_p / c_v`` in three space dimensions.
Here, ``\rho`` is the density, ``\vec{v} = (v_1, v_2, v_3)`` the velocity, ``\theta`` the
potential temperature, ``\phi`` the gravitational potential,
```math
\pi = \left( \frac{p}{p_0} \right)^{R / c_p} = \left( \frac{R \rho \theta}{p_0} \right)^{R / c_v}
```
the Exner pressure, and
```math
p = p_0 \left( \frac{R \rho \theta}{p_0} \right)^\gamma
```
the pressure, where ``p_0 = 10^5\,\mathrm{Pa}`` is the reference pressure and
``R = c_p - c_v`` the gas constant.
"""
struct CompressibleEulerVectorInvariantEquationsWithGravity3D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{3, 6}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
    function CompressibleEulerVectorInvariantEquationsWithGravity3D(; c_p, c_v, gravity)
        c_p, c_v, g = promote(c_p, c_v, gravity)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        K = p_0 * (R / p_0)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        return new{typeof(c_p)}(p_0, c_p, c_v, g, R,
                                gamma,
                                inv_gamma_minus_one,
                                K, stolarsky_factor)
    end
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    ("rho", "v1", "v2", "v3", "rho_theta", "phi")
end

varnames(::typeof(cons2prim), ::CompressibleEulerVectorInvariantEquationsWithGravity3D) = ("rho",
                                                                                           "v1",
                                                                                           "v2",
                                                                                           "v3",
                                                                                           "p",
                                                                                           "phi")

have_nonconservative_terms(::CompressibleEulerVectorInvariantEquationsWithGravity3D) = Trixi.True()

# The numerical fluxes of these equations combine the conservative and the nonconservative
# contributions. Thus, the boundary condition returns only the first contribution of the
# surface flux.
@inline function boundary_condition_slip_wall(u_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function::Union{typeof(flux_etec),
                                                                           typeof(flux_lmars_es)},
                                              equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4] - 2 * u_normal * normal[3],
                         u_inner[5], u_inner[6])

    # calculate the boundary flux
    flux, _ = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

    return flux
end

"""
    flux_etec(u_ll, u_rr, normal_direction,
              equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)

Entropy and total energy conservative two-point flux by
- Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
  Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
  Euler Equations in Vector-Invariant Form
  [arXiv: 2610.05820](https://arxiv.org/abs/2610.05820)

It returns the conservative and the nonconservative contribution of the flux as a tuple,
see [`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref).
"""
@inline function flux_etec(u_ll, u_rr, normal_direction::AbstractVector,
                           equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, rho_theta_ll, phi_ll = u_ll
    rho_rr, v1_rr, v2_rr, v3_rr, rho_theta_rr, phi_rr = u_rr
    exner_ll = exner_pressure(u_ll, equations)
    exner_rr = exner_pressure(u_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]

    # Mean values
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kinetic_energy_avg = 0.25f0 * (v1_ll^2 + v2_ll^2 + v3_ll^2 +
                          v1_rr^2 + v2_rr^2 + v3_rr^2)
    v1_avg = 0.5f0 * (rho_ll * v1_ll + rho_rr * v1_rr) / rho_avg
    v2_avg = 0.5f0 * (rho_ll * v2_ll + rho_rr * v2_rr) / rho_avg
    v3_avg = 0.5f0 * (rho_ll * v3_ll + rho_rr * v3_rr) / rho_avg
    theta_mean = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)

    # Conservative part
    f1 = 0.5f0 * (rho_ll * v_dot_n_ll + rho_rr * v_dot_n_rr)
    f2 = kinetic_energy_avg * normal_direction[1]
    f3 = kinetic_energy_avg * normal_direction[2]
    f4 = kinetic_energy_avg * normal_direction[3]
    f5 = f1 * theta_mean

    # Nonconservative part, i.e., the vorticity, pressure gradient, and gravity terms.
    # Since it is antisymmetric, the right contribution is obtained by flipping the sign.
    gradient = equations.c_p * theta_mean * (exner_rr - exner_ll) + (phi_rr - phi_ll)
    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    vorticity_1 = v2_avg *
                  (jump_v1 * normal_direction[2] - jump_v2 * normal_direction[1]) +
                  v3_avg *
                  (jump_v1 * normal_direction[3] - jump_v3 * normal_direction[1])
    vorticity_2 = v1_avg *
                  (jump_v2 * normal_direction[1] - jump_v1 * normal_direction[2]) +
                  v3_avg *
                  (jump_v2 * normal_direction[3] - jump_v3 * normal_direction[2])
    vorticity_3 = v1_avg *
                  (jump_v3 * normal_direction[1] - jump_v1 * normal_direction[3]) +
                  v2_avg *
                  (jump_v3 * normal_direction[2] - jump_v2 * normal_direction[3])
    g2 = vorticity_1 + gradient * normal_direction[1]
    g3 = vorticity_2 + gradient * normal_direction[2]
    g4 = vorticity_3 + gradient * normal_direction[3]

    zero_ = zero(eltype(u_ll))
    flux_left = SVector(f1, f2 + 0.5f0 * g2, f3 + 0.5f0 * g3, f4 + 0.5f0 * g4, f5,
                        zero_)
    flux_right = SVector(f1, f2 - 0.5f0 * g2, f3 - 0.5f0 * g3, f4 - 0.5f0 * g4, f5,
                         zero_)

    return flux_left, flux_right
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_etec), ::CompressibleEulerVectorInvariantEquationsWithGravity3D) = Trixi.True()

"""
    flux_lmars_es(u_ll, u_rr, normal_direction,
                  equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)

Entropy-stable and well-balanced version of the low Mach number approximate Riemann solver
[`FluxLMARS`](@ref) by
- Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
  Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
  Euler Equations in Vector-Invariant Form
  [arXiv: 2610.05820](https://arxiv.org/abs/2610.05820)

Like [`flux_etec`](@ref), it returns the conservative and the nonconservative contribution
of the surface flux as a tuple, see
[`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref).

The geopotential is continuous across interfaces, so the nonconservative gravity term, which
is proportional to the jump of the geopotential, vanishes on the surface and is not computed
here.
"""
@inline function flux_lmars_es(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, rho_theta_ll = u_ll
    rho_rr, v1_rr, v2_rr, v3_rr, rho_theta_rr = u_rr
    exner_ll = exner_pressure(u_ll, equations)
    exner_rr = exner_pressure(u_rr, equations)
    theta_ll = rho_theta_ll / rho_ll
    theta_rr = rho_theta_rr / rho_rr

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]

    norm_ = norm(normal_direction)
    inv_norm = inv(norm_)

    # Arithmetic mean of the local speeds of sound, where T = θ π
    c_ll = sqrt(equations.gamma * equations.R * theta_ll * exner_ll)
    c_rr = sqrt(equations.gamma * equations.R * theta_rr * exner_rr)
    c_avg = 0.5f0 * (c_ll + c_rr)

    # Mean values
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kinetic_energy_avg = 0.25f0 * (v1_ll^2 + v2_ll^2 + v3_ll^2 +
                          v1_rr^2 + v2_rr^2 + v3_rr^2)
    v1_avg = 0.5f0 * (rho_ll * v1_ll + rho_rr * v1_rr) / rho_avg
    v2_avg = 0.5f0 * (rho_ll * v2_ll + rho_rr * v2_rr) / rho_avg
    v3_avg = 0.5f0 * (rho_ll * v3_ll + rho_rr * v3_rr) / rho_avg
    theta_mean = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    theta_grad_exner = equations.c_p * theta_mean * (exner_rr - exner_ll)

    # Upwinded mass flux and entropy-stable flux of the potential temperature: the
    # entropy-conservative mean `1 / {{1/theta}}_ln` of `flux_etec` and an upwind correction
    # scaled by the jump
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -
                  0.5f0 / c_avg * theta_grad_exner * norm_
    if v_interface > 0
        f1 = rho_ll * v_interface
        sign_v = one(f1)
    else
        f1 = rho_rr * v_interface
        sign_v = -one(f1)
    end
    theta_min = min(theta_ll, theta_rr)
    f5 = f1 *
         (theta_mean - theta_min / (2 * theta_avg) * (theta_rr - theta_ll) * sign_v)

    # Momentum fluxes with dissipation of the normal velocity scaled by the speed of sound
    # and of all velocity components scaled by the advection speed
    diss_normal = 0.5f0 * c_avg / rho_avg *
                  (rho_rr * v_dot_n_rr - rho_ll * v_dot_n_ll) * inv_norm
    diss_advection = 0.25f0 * abs(v_dot_n_ll + v_dot_n_rr) / rho_avg
    f2 = (kinetic_energy_avg - diss_normal) * normal_direction[1] -
         diss_advection * (rho_rr * v1_rr - rho_ll * v1_ll)
    f3 = (kinetic_energy_avg - diss_normal) * normal_direction[2] -
         diss_advection * (rho_rr * v2_rr - rho_ll * v2_ll)
    f4 = (kinetic_energy_avg - diss_normal) * normal_direction[3] -
         diss_advection * (rho_rr * v3_rr - rho_ll * v3_ll)

    # Nonconservative part, i.e., the vorticity and pressure gradient terms
    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    vorticity_1 = v2_avg *
                  (jump_v1 * normal_direction[2] - jump_v2 * normal_direction[1]) +
                  v3_avg *
                  (jump_v1 * normal_direction[3] - jump_v3 * normal_direction[1])
    vorticity_2 = v1_avg *
                  (jump_v2 * normal_direction[1] - jump_v1 * normal_direction[2]) +
                  v3_avg *
                  (jump_v2 * normal_direction[3] - jump_v3 * normal_direction[2])
    vorticity_3 = v1_avg *
                  (jump_v3 * normal_direction[1] - jump_v1 * normal_direction[3]) +
                  v2_avg *
                  (jump_v3 * normal_direction[2] - jump_v2 * normal_direction[3])
    g2 = vorticity_1 + theta_grad_exner * normal_direction[1]
    g3 = vorticity_2 + theta_grad_exner * normal_direction[2]
    g4 = vorticity_3 + theta_grad_exner * normal_direction[3]

    zero_ = zero(eltype(u_ll))
    flux_left = SVector(f1, f2 + 0.5f0 * g2, f3 + 0.5f0 * g3, f4 + 0.5f0 * g4, f5,
                        zero_)
    flux_right = SVector(f1, f2 - 0.5f0 * g2, f3 - 0.5f0 * g3, f4 - 0.5f0 * g4, f5,
                         zero_)

    return flux_left, flux_right
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_lmars_es), ::CompressibleEulerVectorInvariantEquationsWithGravity3D) = Trixi.True()

@inline function Trixi.prim2cons(prim,
                                 equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho, v1, v2, v3, p, phi = prim
    rho_theta = (p / equations.K)^(1 / equations.gamma)
    return SVector(rho, v1, v2, v3, rho_theta, phi)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho, v1, v2, v3, rho_theta, phi = u
    p = equations.K * rho_theta^equations.gamma
    return SVector(rho, v1, v2, v3, p, phi)
end

# Exner pressure π = (p / p_0)^(R / c_p) = (R ρθ / p_0)^(R / c_v)
@inline function exner_pressure(u,
                                equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho_theta = u[5]
    return (equations.R * rho_theta / equations.p_0)^(equations.R / equations.c_v)
end

# Entropy variables w = ∂U/∂u of the entropy U = -ρ s / (γ - 1) with s = log(p) - γ log(ρ),
# see `entropy`. They vanish for the velocity, since U does not depend on it.
@inline function cons2entropy(u,
                              equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho, _, _, _, rho_theta = u

    w1 = -(log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma) *
         equations.inv_gamma_minus_one
    w5 = -equations.gamma * rho / rho_theta * equations.inv_gamma_minus_one

    zero_ = zero(eltype(u))
    return SVector(w1, zero_, zero_, zero_, w5, zero_)
end

@inline function entropy(cons,
                         equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    p = equations.K * cons[5]^equations.gamma
    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])
    S = -s * cons[1] / (equations.gamma - 1)
    return S
end

@inline function pressure(cons,
                          equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    p = equations.K * exp(equations.gamma * log(cons[5]))
    return p
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    rho, v1, v2, v3, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c, abs(v3) + c
end
end # @muladd

# Specializations of the `FluxTurbo` interface of Trixi.jl.
# We precompute `rho, v1, v2, v3, theta, exner, phi, log(theta)`, i.e., two quantities more
# than the number of variables.
@inline function Trixi.nturbovars(::typeof(flux_etec),
                                  ::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    return Val(8)
end

@inline function Trixi.cons2turbo(::typeof(flux_etec),
                                  rho, v1, v2, v3, rho_theta, phi,
                                  equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    theta = rho_theta / rho
    exner = exp(equations.R / equations.c_v *
                log(equations.R * rho_theta / equations.p_0))

    return rho, v1, v2, v3, theta, exner, phi, log(theta)
end

@inline function Trixi.flux_turbo(::typeof(flux_etec),
                                  rho_ll, v1_ll, v2_ll, v3_ll, theta_ll, exner_ll, phi_ll,
                                  log_theta_ll,
                                  rho_rr, v1_rr, v2_rr, v3_rr, theta_rr, exner_rr, phi_rr,
                                  log_theta_rr,
                                  normal_direction_1, normal_direction_2,
                                  normal_direction_3,
                                  equations::CompressibleEulerVectorInvariantEquationsWithGravity3D)
    v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
                 v3_ll * normal_direction_3
    v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
                 v3_rr * normal_direction_3

    # Mean values
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kinetic_energy_avg = 0.25f0 * (v1_ll^2 + v2_ll^2 + v3_ll^2 +
                          v1_rr^2 + v2_rr^2 + v3_rr^2)
    v1_avg = 0.5f0 * (rho_ll * v1_ll + rho_rr * v1_rr) / rho_avg
    v2_avg = 0.5f0 * (rho_ll * v2_ll + rho_rr * v2_rr) / rho_avg
    v3_avg = 0.5f0 * (rho_ll * v3_ll + rho_rr * v3_rr) / rho_avg

    # We inline the logarithmic mean to allow LoopVectorization.jl to optimize
    # it efficiently. This is equivalent to
    # theta_mean = inv_ln_mean(1 / theta_ll, 1 / theta_rr)
    #            = theta_ll * theta_rr * inv_ln_mean(theta_ll, theta_rr)
    RealT = typeof(rho_ll)
    epsilon_z = convert(RealT, 1.0e-4)
    # Coefficients of the Taylor series expansion used for the logarithmic mean
    c1 = convert(RealT, 2 / 3)
    c2 = convert(RealT, 2 / 5)
    c3 = convert(RealT, 2 / 7)
    theta_sum = theta_ll + theta_rr
    theta_jump = theta_rr - theta_ll
    z = theta_jump^2 / theta_sum^2
    special_path = (2 + z * (c1 + z * (c2 + c3 * z))) / theta_sum
    regular_path = (log_theta_rr - log_theta_ll) / theta_jump
    theta_mean = theta_ll * theta_rr * ifelse(z < epsilon_z, special_path, regular_path)

    # Conservative part, i.e., `flux_etec`
    f1 = 0.5f0 * (rho_ll * v_dot_n_ll + rho_rr * v_dot_n_rr)
    f2 = kinetic_energy_avg * normal_direction_1
    f3 = kinetic_energy_avg * normal_direction_2
    f4 = kinetic_energy_avg * normal_direction_3
    f5 = f1 * theta_mean

    # Nonconservative part. Since it is antisymmetric, the right contribution is obtained
    # by flipping the sign.
    gradient = equations.c_p * theta_mean * (exner_rr - exner_ll) + (phi_rr - phi_ll)
    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    vorticity_1 = v2_avg * (jump_v1 * normal_direction_2 - jump_v2 * normal_direction_1) +
                  v3_avg * (jump_v1 * normal_direction_3 - jump_v3 * normal_direction_1)
    vorticity_2 = v1_avg * (jump_v2 * normal_direction_1 - jump_v1 * normal_direction_2) +
                  v3_avg * (jump_v2 * normal_direction_3 - jump_v3 * normal_direction_2)
    vorticity_3 = v1_avg * (jump_v3 * normal_direction_1 - jump_v1 * normal_direction_3) +
                  v2_avg * (jump_v3 * normal_direction_2 - jump_v2 * normal_direction_3)
    g2 = vorticity_1 + gradient * normal_direction_1
    g3 = vorticity_2 + gradient * normal_direction_2
    g4 = vorticity_3 + gradient * normal_direction_3

    zero_ = zero(RealT)
    flux_left = (f1, f2 + 0.5f0 * g2, f3 + 0.5f0 * g3, f4 + 0.5f0 * g4, f5, zero_)
    flux_right = (f1, f2 - 0.5f0 * g2, f3 - 0.5f0 * g3, f4 - 0.5f0 * g4, f5, zero_)

    return flux_left, flux_right
end
