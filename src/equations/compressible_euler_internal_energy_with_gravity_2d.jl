# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw""" 
CompressibleEulerInternalEnergyEquationsWithGravity2D(; c_p, c_v, gravity) 

The compressible Euler equations with gravity
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho v_2 \\ \rho e_{\text{internal}}
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
 \rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ (\rho e_{\text{internal}} +p) v_1
\end{pmatrix}
+
\frac{\partial}{\partial y}
\begin{pmatrix}
\rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ (\rho e_{\text{internal}} +p) v_2
\end{pmatrix}
+
\begin{pmatrix}
0 \\ 0 \\ 0 \\ - v_1 \frac{\partial p}{\partial x} - v_2 \frac{\partial p}{\partial y}
\end{pmatrix}
=
\begin{pmatrix}
0 \\ - \rho \frac{\partial}{\partial x} \phi \\ - \rho \frac{\partial}{\partial y} \phi \\
\end{pmatrix}
```
for an ideal gas with ratio of specific heats gamma in two space dimensions. Here, ``\rho`` is the density, ``v_1, v_2`` are the velocities, ``e_{\text{internal}}`` is the specific internal energy, ``\phi`` is the gravitational potential, and
```math
p = (\gamma - 1) \rho e_{\text{internal}}
```
the pressure.
"""
struct CompressibleEulerInternalEnergyEquationsWithGravity2D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{2, 5}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    function CompressibleEulerInternalEnergyEquationsWithGravity2D(; c_p, c_v,
                                                                   gravity)
        c_p, c_v, g = promote(c_p, c_v, gravity)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        return new{typeof(c_p)}(p_0, c_p, c_v, g, R,
                                gamma,
                                inv_gamma_minus_one)
    end
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    ("rho", "rho_v1", "rho_v2", "rho_e_internal", "phi")
end

varnames(::typeof(cons2prim),
::CompressibleEulerInternalEnergyEquationsWithGravity2D) = ("rho",
                                                            "v1",
                                                            "v2",
                                                            "p", "phi")

have_nonconservative_terms(::CompressibleEulerInternalEnergyEquationsWithGravity2D) = True()

@inline function boundary_condition_slip_wall(u_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_functions,
                                              equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions
    # compute the normal momentum
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4], u_inner[5])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

"""
	flux_conservative_etec(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)

Entropy conserving, total energy conserving and kinetic energy preserving two-point flux by
-  Marco Artiano, Hendrik Ranocha (2026)
   On Affordable High-Order Entropy-Conservative/Stable and 
   Well-Balanced Methods for Nonconservative Hyperbolic Systems
   [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
@inline function flux_conservative_etec(u_ll, u_rr,
                                        normal_direction::AbstractVector,
                                        equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)

    p_avg = 0.5f0 * (p_ll + p_rr)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * inv_rho_p_mean * equations.inv_gamma_minus_one + p_avg * v_dot_n_avg
    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

"""
	flux_nonconservative_etec(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)

Nonconservative part of the entropy conserving, total energy conserving and kinetic energy preserving two-point flux by
-  Marco Artiano, Hendrik Ranocha (2026)
   On Affordable High-Order Entropy-Conservative/Stable and 
   Well-Balanced Methods for Nonconservative Hyperbolic Systems
  [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
@inline function flux_nonconservative_etec(u_ll, u_rr,
                                           normal_direction::AbstractVector,
                                           equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    v_dot_n_avg = 0.5f0 * (v_dot_n_rr + v_dot_n_ll)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)

    phi_jump = phi_rr - phi_ll
    gravity = rho_mean * phi_jump
    v_p_mean = -v_dot_n_avg * (p_rr - p_ll)
    return SVector(zero(eltype(u_ll)), gravity * normal_direction[1],
                   gravity * normal_direction[2], v_p_mean, zero(eltype(u_ll)))
end

"""
	flux_conservative_es(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)

Entropy stable two-point flux by
-  Marco Artiano, Hendrik Ranocha (2026)
   On Affordable High-Order Entropy-Conservative/Stable and 
   Well-Balanced Methods for Nonconservative Hyperbolic Systems
   [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
@inline function flux_conservative_es(u_ll, u_rr, normal_direction::AbstractVector,
                                      equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    v_dot_n_avg = 0.5f0 * (v_dot_n_rr + v_dot_n_ll)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)
    a = 0.5f0 * (c_ll + c_rr)
    norm_ = norm(normal_direction)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    p_interface = 0.5f0 * (p_ll + p_rr) -
                  0.5f0 * a * rho_avg * (v_dot_n_rr - v_dot_n_ll) / norm_
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -
                  1 / (2 * a * rho_avg) * (p_rr - p_ll) * norm_
    p_avg = 0.5f0 * (p_ll + p_rr)
    if (v_interface >= 0)
        rho_upwind = rho_mean - 0.5f0 * (rho_rr - rho_ll)
        _, f2, f3, _, _ = u_ll * v_interface
    else
        rho_upwind = rho_mean + 0.5f0 * (rho_rr - rho_ll)
        _, f2, f3, _, _ = u_rr * v_interface
    end
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)

    # Calculate fluxes depending on normal_direction
    f1 = rho_upwind * v_interface
    f2 = f2 + p_interface * normal_direction[1]
    f3 = f3 + p_interface * normal_direction[2]
    f4 = f1 * inv_rho_p_mean * equations.inv_gamma_minus_one + p_avg * v_interface

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

"""
	flux_nonconservative_es(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)

Nonconservative part of the entropy stable two-point flux by
-  Marco Artiano, Hendrik Ranocha (2026)
   On Affordable High-Order Entropy-Conservative/Stable and 
   Well-Balanced Methods for Nonconservative Hyperbolic Systems
   [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
@inline function flux_nonconservative_es(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    v_dot_n_avg = 0.5f0 * (v_dot_n_rr + v_dot_n_ll)

    # Compute the necessary mean values
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)
    a = 0.5f0 * (c_ll + c_rr)
    norm_ = norm(normal_direction)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -
                  1 / (2 * a * rho_avg) * (p_rr - p_ll) * norm_
    p_avg = 0.5f0 * (p_ll + p_rr)
    v_p_mean = -v_interface * (p_rr - p_ll)

    return SVector(zero(eltype(u_ll)), zero(eltype(u_ll)), zero(eltype(u_ll)), v_p_mean,
                   zero(eltype(u_ll)))
end

@inline function prim2cons(prim,
                           equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, v1, v2, p, phi = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_e_internal = p * equations.inv_gamma_minus_one
    return SVector(rho, rho_v1, rho_v2, rho_e_internal, phi)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e_internal, phi = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * rho_e_internal
    return SVector(rho, v1, v2, p, phi)
end

@inline function cons2cons(u,
                           equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    return u
end

@inline function cons2entropy(u,
                              equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e_internal, phi = u

    w1 = log(rho_e_internal * (equations.gamma - 1) / rho^equations.gamma) -
         equations.gamma
    w4 = rho / (rho_e_internal * (equations.gamma - 1))

    return SVector(w1, zero(eltype(u)), zero(eltype(u)), w4, zero(eltype(u)))
end

@inline function entropy(cons,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    p = pressure(cons, equations)
    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])
    S = -s * cons[1] / (equations.gamma - 1)
    return S
end

@inline function pressure(cons,
                          equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    p = (equations.gamma - 1) * cons[4]
    return p
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll, phi = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speeds
    # left
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_, abs(v_rr) + c_rr * norm_)
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, v1, v2, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c
end
end # @muladd
