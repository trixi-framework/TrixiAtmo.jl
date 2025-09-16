using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, FluxLMARS,
              boundary_condition_slip_wall, pressure, energy_total, energy_kinetic

@muladd begin
#! format: noindent
struct CompressibleEulerPotentialTemperatureEquations3D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{3, 5}
    p_0::RealT
    c_p::RealT
    c_v::RealT
    g::RealT
    R::RealT
    gamma::RealT
    inv_gamma_minus_one::RealT
    K::RealT
    stolarsky_factor::RealT
end

function CompressibleEulerPotentialTemperatureEquations3D(; g = 9.81, RealT = Float64)
    p_0 = 100_000.0
    c_p = 1004.0
    c_v = 717.0
    R = c_p - c_v
    gamma = c_p / c_v
    inv_gamma_minus_one = inv(gamma - 1)
    K = p_0 * (R / p_0)^gamma
    stolarsky_factor = (gamma - 1.0) / gamma
    return CompressibleEulerPotentialTemperatureEquations3D{RealT}(p_0, c_p, c_v, g, R,
                                                                   gamma,
                                                                   inv_gamma_minus_one,
                                                                   K, stolarsky_factor)
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerPotentialTemperatureEquations3D)
    ("rho", "rho_v1", "rho_v2", "rho_v3", "rho_theta")
end

varnames(::typeof(cons2prim),
         ::CompressibleEulerPotentialTemperatureEquations3D) = ("rho",
                                                                "v1",
                                                                "v2",
                                                                "v3",
                                                                "p1")

# Calculate 1D flux for a single point in the normal direction.
# Note, this directional vector is not normalized.
@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
               v3 * normal_direction[3]
    rho_v_normal = rho * v_normal
    f1 = rho_v_normal
    f2 = rho_v_normal * v1 + p * normal_direction[1]
    f3 = rho_v_normal * v2 + p * normal_direction[2]
    f4 = rho_v_normal * v3 + p * normal_direction[3]
    f5 = rho_theta * v_normal
    return SVector(f1, f2, f3, f4, f5)
end

@inline function flux(u, orientation::Integer,
                      equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = equations.K * rho_theta^equations.gamma
    if orientation == 1
        f1 = rho_v1
        f2 = rho_v1 * v1 + p
        f3 = rho_v1 * v2
        f4 = rho_v1 * v3
        f5 = rho_theta * v1
    elseif orientation == 2
        f1 = rho_v2
        f2 = rho_v2 * v1
        f3 = rho_v2 * v2 + p
        f4 = rho_v2 * v3
        f5 = rho_theta * v2
    else
        f1 = rho_v3
        f2 = rho_v3 * v1
        f3 = rho_v3 * v2
        f4 = rho_v3 * v3 + p
        f5 = rho_theta * v3
    end

    return SVector(f1, f2, f3, f4, f5)
end

@inline function source_terms_gravity(u, x, t,
                                      equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, _, _, _, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)),
                   zero(eltype(u)), -equations.g * rho, zero(eltype(u)))
end

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerPotentialTemperatureEquations3D)
    a = flux_lmars.speed_of_sound
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f1, f2, f3, f4, f5 = u_ll * v_interface
    else
        f1, f2, f3, f4, f5 = u_rr * v_interface
    end

    return SVector(f1,
                   f2 + p_interface * normal_direction[1],
                   f3 + p_interface * normal_direction[2],
                   f4 + p_interface * normal_direction[3],
                   f5)
end

"""
	flux_tec(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquationsPotentialTemperature1D)

Total energy conservative two-point flux by
-  Artiano et al. (2025), pre-print
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""

@inline function flux_tec(u_ll, u_rr, normal_direction::AbstractVector,
                          equations::CompressibleEulerPotentialTemperatureEquations3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]
    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr
    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)

    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    return SVector(f1, f2, f3, f4, f5)
end

"""
	flux_ec(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquationsPotentialTemperature1D)

Entropy conservative two-point flux by
-  Artiano et al. (2025), pre-print
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""

@inline function flux_ec(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerPotentialTemperatureEquations3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]
    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr
    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    return SVector(f1, f2, f3, f4, f5)
end

"""
	flux_etec(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquationsPotentialTemperature1D)

Entropy and total energy conservative two-point flux by
-  Artiano et al. (2025), pre-print
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""

@inline function flux_etec(u_ll, u_rr, normal_direction::AbstractVector,
                           equations::CompressibleEulerPotentialTemperatureEquations3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]
    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr
    # Compute the necessary mean values

    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f1 = f5 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]

    return SVector(f1, f2, f3, f4, f5)
end

@inline function prim2cons(prim,
                           equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, v1, v2, v3, p = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_v3 = rho * v3
    rho_theta = (p / equations.K)^(1 / equations.gamma)
    return SVector(rho, rho_v1, rho_v2, rho_v3, rho_theta)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = equations.K * rho_theta^equations.gamma
    return SVector(rho, v1, v2, v3, p)
end

@inline function cons2cons(u,
                           equations::CompressibleEulerPotentialTemperatureEquations3D)
    return u
end

@inline function cons2entropy_rhoe(u,
                                   equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u

    w1 = -0.5f0 * (rho_v1^2 + rho_v2^2 + rho_v3^2) / rho^2
    w2 = rho_v1 / rho
    w3 = rho_v2 / rho
    w4 = rho_v3 / rho
    w5 = equations.gamma * equations.inv_gamma_minus_one * equations.K *
         (rho_theta)^(equations.gamma - 1)

    return SVector(w1, w2, w3, w4, w5)
end

@inline function cons2entropy(u,
                              equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u

    w1 = log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma
    w2 = 0.0
    w3 = 0.0
    w4 = 0.0
    w5 = rho / rho_theta * equations.gamma

    return SVector(w1, w2, w3, w4, w5)
end

@inline function energy_total(cons,
                              equations::CompressibleEulerPotentialTemperatureEquations3D)
    p = equations.K * cons[5]^equations.gamma
    U = (p / (equations.gamma - 1) +
         0.5f0 * (cons[2]^2 + cons[3]^2 + cons[4]^2) / (cons[1]))

    return U
end

@inline function entropy(cons,
                         equations::CompressibleEulerPotentialTemperatureEquations3D)
    p = equations.K * cons[5]^equations.gamma
    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])
    S = -s * cons[1] / (equations.gamma - 1.0)
    return S
end

@inline function energy_kinetic(cons,
                                equations::CompressibleEulerPotentialTemperatureEquations3D)
    return 0.5f0 * (cons[2]^2 + cons[3]^2 + cons[4]^2) / (cons[1])
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, v1, v2, v3, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c, abs(v3) + c
end

@inline function density_pressure(u,
                                  equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u
    rho_times_p = rho * equations.p_0 *
                  (equations.R * rho_theta / equations.p_0)^equations.gamma
    return rho_times_p
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    else # orientation == 3
        v_ll = v3_ll
        v_rr = v3_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    λ_max = max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1]
            + v2_ll * normal_direction[2]
            + v3_ll * normal_direction[3])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            + v2_rr * normal_direction[2]
            + v3_rr * normal_direction[3])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

@inline function pressure(cons,
                          equations::CompressibleEulerPotentialTemperatureEquations3D)
    _, _, _, _, p = cons2prim(cons, equations)
    return p
end
end # @muladd
