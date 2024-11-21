# Implemented by Lucas Gemein
# https://github.com/NichtLucas/Trixi.jl/tree/thesis_gemein_2022

using Trixi
using Trixi: ln_mean, inv_ln_mean
import Trixi: varnames, flux_chandrashekar, boundary_condition_slip_wall,
              cons2prim, cons2entropy, max_abs_speed_naive, max_abs_speeds,
              entropy, energy_total, flux, stolarsky_mean

@muladd begin
#! format: noindent
struct CompressibleMoistEulerPotentialTemperatureEquations2D{RealT <: Real} <:
       AbstractCompressibleMoistEulerEquations{2, 6}
    p_0::RealT   # constant reference pressure 1000 hPa(100000 Pa)
    c_pd::RealT   # dry air constant
    c_vd::RealT   # dry air constant
    R_d::RealT   # dry air gas constant
    c_pv::RealT   # moist air constant
    c_vv::RealT   # moist air constant
    R_v::RealT   # moist air gas constant
    c_pl::RealT # liqid water constant
    g::RealT # gravitation constant
    inv_gamma_minus_one::RealT # ratio of the gas constant R_d
    gamma::RealT # = inv(kappa- 1); can be used to write slow divisions as fast multiplications
    L_00::RealT # latent heat of evaporation  at 0 K
    K::RealT
end

function CompressibleMoistEulerPotentialTemperatureEquations2D(; g = 9.81, RealT = Float64)
    p_0 = 100000.0
    c_pd = 1004.0
    c_vd = 717.0
    R_d = c_pd - c_vd
    c_pv = 1885.0
    c_vv = 1424.0
    R_v = c_pv - c_vv
    c_pl = 4186.0
    gamma = c_pd / c_vd # = 1/(1 - kappa)
    inv_gamma_minus_one = inv(1 - gamma)
    K = 
    L_00 = 2.5e6 #+ (c_pl - c_pv)*273.15
    K = p_0 * (R_d / p_0)^gamma
    return CompressibleMoistEulerPotentialTemperatureEquations2D{RealT}(p_0, c_pd, c_vd, R_d, c_pv, c_vv,
                                                    R_v, c_pl, g, inv_gamma_minus_one, gamma, L_00, K)
end

@inline function pressure(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack p_0, R_d, R_v, c_pd, c_pv, c_pl, c_vd, c_vv = equations
    rho, rho_v1, rho_v2, rho_theta, rho_v, rho_l = u

    rho_d = rho - (rho_v + rho_l)

    kappa_M = (R_d * rho_d + R_v * rho_v) / (c_pd * rho_d + c_pv * rho_v + c_pl * rho_l)
    p = p_0 * (R_d * rho_theta / p_0)^(1 / (1 - kappa_M))

    return p
end

# Calculate 1D flux for a single point.
@inline function flux(u, orientation::Integer,
                      equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_qc = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    qv = rho_qv / rho
    qc = rho_qc / rho
    p = pressure(u, equations)
    if orientation == 1
        f1 = rho_v1
        f2 = rho_v1 * v1 + p
        f3 = rho_v1 * v2
        f4 = rho_theta * v1
        f5 = rho_v1 * qv
        f6 = rho_v1 * qc
    else
        f1 = rho_v2
        f2 = rho_v2 * v1
        f3 = rho_v2 * v2 + p
        f4 = rho_theta * v2
        f5 = rho_v2 * qv
        f6 = rho_v2 * qc
    end
    return SVector(f1, f2, f3, f4, f5, f6)
end

# Calculate 1D flux for a single point in the normal direction.
# Note, this directional vector is not normalized.
@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_qc = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    qv = rho_qv / rho
    qc = rho_qc / rho
    p = pressure(u, equations)
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    rho_v_normal = rho * v_normal
    f1 = rho_v_normal
    f2 = (rho_v_normal) * v1 + p * normal_direction[1]
    f3 = (rho_v_normal) * v2 + p * normal_direction[2]
    f4 = (rho_theta) * v_normal
    f5 = rho_v_normal * qv
    f6 = (rho_v_normal) * qc
    return SVector(f1, f2, f3, f4, f5, f6)
end

# Slip-wall boundary condition
# Determine the boundary numerical surface flux for a slip wall condition.
# Imposes a zero normal velocity at the wall.
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function,
                                              equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack c_pd, c_pv, c_pl, c_vd, c_vv = equations
    norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    # rotate the internal solution state
    u_local = rotate_to_x(u_inner, normal, equations)

    # compute the primitive variables
    rho_local, v_normal, v_tangent, p_local, qv_local, ql_local = cons2prim(u_local,
                                                                            equations)
    qd_local = 1 - qv_local - ql_local
    gamma = (qd_local * c_pd + qv_local * c_pv + ql_local * c_pl) *
            inv(qd_local * c_vd + qv_local * c_vv + ql_local * c_pl)
    # Get the solution of the pressure Riemann problem
    # See Section 6.3.3 of
    # Eleuterio F. Toro (2009)
    # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
    # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
    if v_normal <= 0.0
        sound_speed = sqrt(gamma * p_local / rho_local) # local sound speed
        p_star = p_local *
                 (1.0 + 0.5 * (gamma - 1) * v_normal / sound_speed)^(2.0 * gamma *
                                                                     inv(gamma - 1))
    else # v_normal > 0.0
        A = 2.0 / ((gamma + 1) * rho_local)
        B = p_local * (gamma - 1) / (gamma + 1)
        p_star = p_local +
                 0.5 * v_normal / A *
                 (v_normal + sqrt(v_normal^2 + 4.0 * A * (p_local + B)))
    end

    # For the slip wall we directly set the flux as the normal velocity is zero
    return SVector(zero(eltype(u_inner)),
                   p_star * normal[1],
                   p_star * normal[2],
                   zero(eltype(u_inner)),
                   zero(eltype(u_inner)),
                   zero(eltype(u_inner))) * norm_
end

# Fix sign for structured mesh.
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                              direction, x, t,
                                              surface_flux_function,
                                              equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    # flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
    # to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
    if isodd(direction)
        boundary_flux = -boundary_condition_slip_wall(u_inner, -normal_direction,
                                                      x, t, surface_flux_function,
                                                      equations)
    else
        boundary_flux = boundary_condition_slip_wall(u_inner, normal_direction,
                                                     x, t, surface_flux_function,
                                                     equations)
    end

    return boundary_flux
end

# Rotate momentum flux. The same as in compressible Euler.
@inline function rotate_to_x(u, normal_vector::AbstractVector,
                             equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    # Apply the 2D rotation matrix with normal and tangent directions of the form
    # [ 1    0    0   0;
    #   0   n_1  n_2  0;
    #   0   t_1  t_2  0;
    #   0    0    0   1 ]
    # where t_1 = -n_2 and t_2 = n_1

    return SVector(u[1],
                   c * u[2] + s * u[3],
                   -s * u[2] + c * u[3],
                   u[4], u[5], u[6])
end

# Recreates the convergence test initial condition from compressible euler 2D.
function initial_condition_warm_bubble(x, t,
                                                equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack R_d, R_v, c_vd, c_vv, c_pl, L_00 = equations
    xc = 0
    zc = 2000
    r = sqrt((x[1] - xc)^2 + (x[2] - zc)^2)
    rc = 2000
    theta_ref = 300
    qv_ref = 0
    qc_ref = 0
    Δtheta = 0
    Δqv = 0
    Δqc = 0
  
    if r <= rc
       Δtheta = 2 * cospi(0.5*r/rc)^2
    end
  
    #Perturbed state:
    theta = theta_ref + Δtheta # potential temperature
    π_exner = 1 - equation._grav / (equation.c_pd * theta) * x[2] # exner pressure
    rho = equation.p_0 / (equation.R_d * theta) * (π_exner)^(equation.c_vd / equation.R_d) # density
  
    v1 = 20
    v2 = 0
    qv = 0
    qc = 0
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_theta = rho * theta
    rho_qv = rho * qv
    rho_qc = rho * qc

    return SVector(rho, rho_v1, rho_v2, rho_e, rho_qv, rho_qc)
end



# Gravity source term
@inline function source_terms_gravity(u, x, t,
                                           equations::CompressibleMoistEulerEquations2D)
    @unpack g = equations
    rho, _, _, _, _, _ = u

    return SVector(zero(eltype(u)), zero(eltype(u)),
                   -g * rho,
                   zero(eltype(u)), zero(eltype(u)), zero(eltype(u)))
end

# Rayleigh damping sponge source term form A. Sridhar et al.,
# Large-eddy simulations with ClimateMachine: a new open-sourcecode for
# atmospheric simulations on GPUs and CPUs, 2 Oct 2021, doi: 10.5194/gmd-15-6259-2022,
# https://arxiv.org/abs/2110.00853 [physics.ao-ph] .
@inline function source_terms_nonhydrostatic_rayleigh_sponge(u, x, t,
                                                             equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_ql = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    z = x[2]

    # relaxed background velocity
    vr1, vr2 = (10.0, 0.0)
    # damping threshold
    z_s = 9000.0
    # boundary top
    z_top = 16000.0
    # positive even power with default value 2
    gamma = 2.0
    # relaxation coefficient > 0
    alpha = 0.5

    tau_s = zero(eltype(u))
    if z > z_s
        tau_s = alpha * sin(0.5 * (z - z_s) * inv(z_top - z_s))^(gamma)
    end

    return SVector(zero(eltype(u)),
                   -tau_s * rho * (v1 - vr1),
                   -tau_s * rho * (v2 - vr2),
                   zero(eltype(u)), zero(eltype(u)), zero(eltype(u)))
end

function source_terms_moist_bubble(u, x, t, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)

    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_qc = u
    RelCloud = 1
    rho_d = rho - rho_qv - rho_qc
    c_pml = equations.c_pd * rho_d + equations.c_pv * rho_qv + equations.c_pl * rho_qc
    c_vml = equations.c_vd * rho_d + equations.c_vv * rho_qv + equations.c_pl * rho_qc
    R_m   = equations.R_d * rho_d + equations.R_v * rho_qv
    kappa_M = R_m / c_pml
    p = pressure(u, equations)
    T = p / R_m
    T_C = T - 273.15
    p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
    a = p_vs / (equations.R_v * T) - rho_qv
    b = rho_qc
    rho_q_cond = RelCloud * (a + b - sqrt(a * a + b * b))
    L = equations.L_00 - (equations.c_pl - equations.c_pv) * T
    du4 = rho_theta * ((-L / (c_pml * T) - log(p / equations.p_0) * kappa_M * (equations.R_v / R_m - equations.c_pv / c_pml) + 
    equations.R_v / R_m) * rho_q_cond +  (log(p / equations.p_0) * kappa_M * (equations.c_pl / c_pml)) * (-rho_q_cond))

    du5 =  rho_q_cond
    du6 = -rho_q_cond

    return SVector(zero(eltype(u)), zero(eltype(u)), -equations.g*rho, du4, du5, du6)
end


# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.
@inline function flux_LMARS(u_ll, u_rr, normal_direction::AbstractVector,
                            equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    a = 360.0
    # Unpack left and right state
    rho_ll, rho_v1_ll, rho_v2_ll, rho_theta_ll, rho_qv_ll, rho_ql_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_theta_rr, rho_qv_rr, rho_ql_rr = u_rr
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)
    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v1_rr = rho_v1_rr / rho_rr
    v2_rr = rho_v2_rr / rho_rr

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # diffusion parameter <= 1
    beta = 1

    # Compute the necessary interface flux components
    norm_ = norm(normal_direction)

    rho = 0.5 * (rho_ll + rho_rr)
    p_interface = 0.5 * (p_ll + p_rr) - beta * 0.5 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5 * (v_ll + v_rr) - beta * 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f1, f2, f3, f4, f5, f6 = u_ll * v_interface
    else
        f1, f2, f3, f4, f5, f6 = u_rr * v_interface
    end

    return SVector(f1,
                   f2 + p_interface * normal_direction[1],
                   f3 + p_interface * normal_direction[2],
                   f4, f5, f6)
end

# Convert conservative variables to primitive.
@inline function cons2prim(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_ql = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p =  pressure(u, equations)
    qv = rho_qv / rho
    ql = rho_ql / rho

    return SVector(rho, v1, v2, p, qv, ql)
end

# Convert conservative variables to entropy
@inline function cons2entropy(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack R_d, R_v, c_pd, c_pv, c_pl, L_00 = equations
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_qc = u

  v1 = rho_v1 / rho
  v2 = rho_v2 / rho
  v_square = v1^2 + v2^2
  p = pressure(u, equations)
  s = log(p) - equations.gamma*log(rho)
  rho_p = rho / p

  w1 = (equations.gamma - s) / (equations.gamma-1) - 0.5 * rho_p * v_square
  w2 = rho_p * v1
  w3 = rho_p * v2
  w4 = -rho_p
  w5 = rho_p * rho_qv / rho
  w6 = rho_p * rho_qc / rho

    return SVector(w1, w2, w3, w4, w5, w6)
end

# Convert primitive to conservative variables.
@inline function prim2cons(prim, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho, v1, v2, p, qv, ql = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_qv = rho * qv
    rho_ql = rho * ql
    rho_theta = (p / equations.K)^(1 / equations.gamma) # TODO
    return SVector(rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_ql)
end


@inline function density(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho = u[1]
    return rho
end

@inline function density_dry(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho_qd = u[1] - (u[5] + u[6])
    return rho_qd
end

@inline function density_vapor(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho_qv = u[5]
    return rho_qv
end





# Calculate kinetic energy for a conservative state `cons`.
@inline function energy_kinetic(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_ql = u
    return (rho_v1^2 + rho_v2^2) / (2 * rho)
end


@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack c_pd, c_pv, c_pl, c_vd, c_vv = equations
    rho_ll, v1_ll, v2_ll, p_ll, qv_ll, ql_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, qv_rr, ql_rr = cons2prim(u_rr, equations)
    qd_ll = 1 - qv_ll - ql_ll
    qd_rr = 1 - qv_rr - ql_rr
    # Get the density and gas gamma
    gamma_ll = (qd_ll * c_pd + qv_ll * c_pv + ql_ll * c_pl) *
               inv(qd_ll * c_vd + qv_ll * c_vv + ql_ll * c_pl)
    gamma_rr = (qd_rr * c_pd + qv_rr * c_pv + ql_rr * c_pl) *
               inv(qd_rr * c_vd + qv_rr * c_vv + ql_rr * c_pl)

    # Compute the sound speeds on the left and right
    v_mag_ll = sqrt(v1_ll^2 + v2_ll^2)
    c_ll = sqrt(gamma_ll * p_ll / rho_ll)
    v_mag_rr = sqrt(v1_rr^2 + v2_rr^2)
    c_rr = sqrt(gamma_rr * p_rr / rho_rr)

    return max(v_mag_ll, v_mag_rr) + max(c_ll, c_rr)
end

# Adjusted version of LLF dissipation from compressible euler.
@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack c_pd, c_pv, c_pl, c_vd, c_vv = equations
    rho_ll, v1_ll, v2_ll, p_ll, qv_ll, ql_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, qv_rr, ql_rr = cons2prim(u_rr, equations)
    qd_ll = 1 - qv_ll - ql_ll
    qd_rr = 1 - qv_rr - ql_rr
    # Get the density and gas gamma
    gamma_ll = (qd_ll * c_pd + qv_ll * c_pv + ql_ll * c_pl) *
               inv(qd_ll * c_vd + qv_ll * c_vv + ql_ll * c_pl)
    gamma_rr = (qd_rr * c_pd + qv_rr * c_pv + ql_rr * c_pl) *
               inv(qd_rr * c_vd + qv_rr * c_vv + ql_rr * c_pl)
    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2])
    c_ll = sqrt(gamma_ll * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2])
    c_rr = sqrt(gamma_rr * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Adjusted version of lambda_max from compressible euler.
@inline function max_abs_speeds(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack c_pd, c_pv, c_pl, c_vd, c_vv = equations
    rho, v1, v2, p, qv, ql = cons2prim(u, equations)
    qd = 1 - qv - ql

    gamma = (qd * c_pd + qv * c_pv + ql * c_pl) * inv(qd * c_vd + qv * c_vv + ql * c_pl)
    c = sqrt(gamma * p / rho)

    return (abs(v1) + c, abs(v2) + c)
end

@inline function cons2aeqpot(u, equations::CompressibleMoistEulerPotentialTemperatureEquations2D)
    @unpack c_pd, c_pv, c_pl, R_d, R_v, p_0, L_00 = equations
    rho, rho_v1, rho_v2, rho_theta, rho_qv, rho_ql = u
    rho_d = rho - rho_qv - rho_ql
    p = pressure(u, equations)
    T = p / (equations.R_d * rho_d + equations.R_v * rho_qv)
    p_v = rho_qv * R_v * T
    p_d = p - p_v
    T_C = T - 273.15
    p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
    H = p_v / p_vs
    r_v = rho_qv / rho_d
    r_l = rho_ql / rho_d
    r_t = r_v + r_l
    L_v = L_00 + (c_pv - c_pl) * T
    c_p = c_pd + r_t * c_pl

    # equivalent potential temperature
    aeq_pot = (T * (p_0 / p_d)^(R_d / c_p) * H^(-r_v * R_v / c_p) *
               exp(L_v * r_v * inv(c_p * T)))

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho

    pot1 = rho
    pot2 = v1
    pot3 = v2
    pot4 = aeq_pot
    pot5 = r_v
    pot6 = r_t
    return SVector(pot1, pot2, pot3, pot4, pot5, pot6)
end

varnames(::typeof(cons2aeqpot), ::CompressibleMoistEulerPotentialTemperatureEquations2D) = ("rho", "v1",
                                                                        "v2",
                                                                        "aeqpottemp",
                                                                        "rv", "rt")


varnames(::typeof(cons2cons), ::CompressibleMoistEulerPotentialTemperatureEquations2D) = ("rho", "rho_v1",
                                                                      "rho_v2", "rho_theta",
                                                                      "rho_qv",
                                                                      "rho_ql")
varnames(::typeof(cons2prim), ::CompressibleMoistEulerPotentialTemperatureEquations2D) = ("rho", "v1", "v2",
                                                                      "p", "qv", "ql")
end # @muladd
