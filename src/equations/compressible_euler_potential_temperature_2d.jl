@muladd begin
#! format: noindent
struct CompressibleEulerPotentialTemperatureEquations2D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{2, 4}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg·K)
    R::RealT # gas constant
    gamma::RealT # ratio of specific heats 
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
    function CompressibleEulerPotentialTemperatureEquations2D(; c_p, c_v)
        c_p, c_v = promote(c_p, c_v)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        K = p_0 * (R / p_0)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        return new{typeof(gamma)}(p_0, c_p, c_v, R,
                                  gamma,
                                  inv_gamma_minus_one,
                                  K, stolarsky_factor)
    end
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerPotentialTemperatureEquations2D)
    ("rho", "rho_v1", "rho_v2", "rho_theta")
end

varnames(::typeof(cons2prim),
::CompressibleEulerPotentialTemperatureEquations2D) = ("rho",
                                                       "v1",
                                                       "v2",
                                                       "p1")

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    rho_v_normal = rho * v_normal
    p = pressure(u, equations)

    f1 = rho_v_normal
    f2 = (rho_v_normal) * v1 + p * normal_direction[1]
    f3 = (rho_v_normal) * v2 + p * normal_direction[2]
    f4 = (rho_theta) * v_normal

    return SVector(f1, f2, f3, f4)
end

# Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer,
                      equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = pressure(u, equations)
    if orientation == 1
        f1 = rho_v1
        f2 = rho_v1 * v1 + p
        f3 = rho_v1 * v2
        f4 = rho_theta * v1
    else
        f1 = rho_v2
        f2 = rho_v2 * v1
        f3 = rho_v2 * v2 + p
        f4 = rho_theta * v2
    end
    return SVector(f1, f2, f3, f4)
end

@inline function boundary_condition_slip_wall(u_inner, orientation,
                                              direction, x, t,
                                              surface_flux_function,
                                              equations::CompressibleEulerPotentialTemperatureEquations2D)

    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])
    else # orientation == 2
        u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4])
    end
    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerPotentialTemperatureEquations2D)
    a = flux_lmars.speed_of_sound
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f1, f2, f3, f4 = u_ll * v_interface
    else
        f1, f2, f3, f4 = u_rr * v_interface
    end

    return SVector(f1,
                   f2 + p_interface * normal_direction[1],
                   f3 + p_interface * normal_direction[2],
                   f4)
end

"""
	flux_ec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature2D)

Entropy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_ec(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerPotentialTemperatureEquations2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    _, _, _, rho_theta_ll = u_ll
    _, _, _, rho_theta_rr = u_rr

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)

    return SVector(f1, f2, f3, f4)
end

"""
	flux_tec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature2D)

Total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_tec(u_ll, u_rr, normal_direction::AbstractVector,
                          equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    rho_theta_ll = last(u_ll)
    rho_theta_rr = last(u_rr)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)

    return SVector(f1, f2, f3, f4)
end

"""
	flux_etec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature2D)

Entropy and total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_etec(u_ll, u_rr, normal_direction::AbstractVector,
                           equations::CompressibleEulerPotentialTemperatureEquations2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    _, _, _, rho_theta_ll = u_ll
    _, _, _, rho_theta_rr = u_rr

    # Compute the necessary mean values
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f1 = f4 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]

    return SVector(f1, f2, f3, f4)
end

@inline function prim2cons(prim,
                           equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, v1, v2, p = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
    return SVector(rho, rho_v1, rho_v2, rho_theta)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = equations.K * rho_theta^equations.gamma
    return SVector(rho, v1, v2, p)
end

@inline function cons2cons(u,
                           equations::CompressibleEulerPotentialTemperatureEquations2D)
    return u
end

@inline function cons2entropy(u,
                              equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta = u

    w1 = log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma
    w4 = rho / rho_theta * equations.gamma

    return SVector(w1, zero(eltype(u)), zero(eltype(u)), w4)
end

@inline function entropy_thermodynamic(cons,
                                       equations::CompressibleEulerPotentialTemperatureEquations2D)
    # Mathematical entropy
    p = equations.p_0 * (equations.R * cons[4] / equations.p_0)^equations.gamma

    U = (p / (equations.gamma - 1) + 0.5f0 * (cons[2]^2 + cons[3]^2) / (cons[1]))

    return U
end

# Default entropy is the mathematical entropy
@inline function entropy(cons,
                         equations::CompressibleEulerPotentialTemperatureEquations2D)
    entropy_thermodynamic(cons, equations)
end

@inline function energy_total(cons,
                              equations::CompressibleEulerPotentialTemperatureEquations2D)
    entropy(cons, equations)
end

@inline function energy_kinetic(cons,
                                equations::CompressibleEulerPotentialTemperatureEquations2D)
    return 0.5f0 * (cons[2]^2 + cons[3]^2) / (cons[1])
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, v1, v2, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    λ_max = max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speed
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

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, orientation::Integer,
                               equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

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
    return max(abs(v_ll) + c_ll * norm_,
               abs(v_rr) + c_rr * norm_)
end
@inline function pressure(cons,
                          equations::CompressibleEulerPotentialTemperatureEquations2D)
    p = equations.K * exp(equations.gamma * log(cons[4]))
    return p
end
end # @muladd
