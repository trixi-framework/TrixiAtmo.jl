@muladd begin
#! format: noindent

struct CompressibleEulerPotentialTemperatureEquationsWithGravity1D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{1, 4}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats 
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
end

function CompressibleEulerPotentialTemperatureEquationsWithGravity1D(; g = 9.81,
                                                                     RealT = Float64)
    p_0 = 100_000
    c_p = 1004
    c_v = 717
    R = c_p - c_v
    gamma = c_p / c_v
    inv_gamma_minus_one = inv(gamma - 1)
    K = p_0 * (R / p_0)^gamma
    stolarsky_factor = (gamma - 1) / gamma
    return CompressibleEulerPotentialTemperatureEquationsWithGravity1D{RealT}(p_0, c_p,
                                                                              c_v, g, R,
                                                                              gamma,
                                                                              inv_gamma_minus_one,
                                                                              K,
                                                                              stolarsky_factor)
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    ("rho", "rho_v1", "rho_theta", "phi")
end

varnames(::typeof(cons2prim),
::CompressibleEulerPotentialTemperatureEquationsWithGravity1D) = ("rho", "v1",
                                                                  "p1", "phi")

have_nonconservative_terms(::CompressibleEulerPotentialTemperatureEquationsWithGravity1D) = Trixi.True()

@inline function boundary_condition_slip_wall(u_inner, orientation,
                                              direction, x, t,
                                              surface_flux_functions,
                                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    ## get the appropriate normal vector from the orientation
    u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])
    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
                                                     equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
                                                     equations)
    end

    return flux, noncons_flux
end

@inline function flux(u, orientation::Integer,
                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho, rho_v1, rho_theta = u
    v1 = rho_v1 / rho
    p = equations.K * exp(log(rho_theta^equations.gamma))
    p = pressure(u, equations)
    f1 = rho_v1
    f2 = rho_v1 * v1 + p
    f3 = rho_theta * v1

    return SVector(f1, f2, f3, 0)
end

"""
   flux_nonconservative_waruzewski_etal(u_ll, u_rr, orientation::Integer, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)

Well-balanced gravity term for isothermal background state
-  Maciej Waruszewski and Jeremy E. Kozdon and Lucas C. Wilcox and Thomas H. Gibson and Francis X. Giraldo (2022),
   Entropy stable discontinuous {G}alerkin methods for balance laws 
   in non-conservative form: Applications to the {E}uler equations with gravity
   [DOI: 10.1016/j.jcp.2022.111507](https://doi.org/10.1016/j.jcp.2022.111507)

The well balanced on curvilinear coordinates was proven by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_waruzewski_etal(u_ll, u_rr,
                                                      orientation::Integer,
                                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho_ll, _, _, phi_ll = u_ll
    rho_rr, _, _, phi_rr = u_rr
    rho_avg = ln_mean(rho_ll, rho_rr)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   rho_avg * phi_jump,
                   zero(eltype(u_ll)),
                   zero(eltype(u_ll)))
end

"""
   flux_nonconservative_artiano_etal(u_ll, u_rr, orientation::Integer, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)

Well-balanced gravity term for constant potential temperature background state by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_artiano_etal(u_ll, u_rr,
                                                   orientation::Integer,
                                                   equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho_ll, _, _, phi_ll = u_ll
    rho_rr, _, _, phi_rr = u_rr
    rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   rho_avg * phi_jump,
                   zero(eltype(u_ll)),
                   zero(eltype(u_ll)))
end

"""
   flux_nonconservative_souza_etal(u_ll, u_rr, orientation::Integer, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)

-  Souza et al. 
   The Flux-Differencing Discontinuous {G}alerkin Method Applied to 
   an Idealized Fully Compressible Nonhydrostatic Dry Atmosphere
   [DOI: 10.1029/2022MS003527] (https://doi.org/10.1029/2022MS003527)
"""
@inline function flux_nonconservative_souza_etal(u_ll, u_rr,
                                                 orientation::Integer,
                                                 equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho_ll, _, _, phi_ll = u_ll
    rho_rr, _, _, phi_rr = u_rr
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   rho_avg * phi_jump,
                   zero(eltype(u_ll)),
                   zero(eltype(u_ll)))
end

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
                                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    a = flux_lmars.speed_of_sound
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    rho = 0.5f0 * (rho_ll + rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v1_rr - v1_ll)
    v_interface = 0.5f0 * (v1_ll + v1_rr) - 1 / (2 * a * rho) * (p_rr - p_ll)

    if (v_interface > 0)
        f1, f2, f3 = u_ll * v_interface
    else
        f1, f2, f3 = u_rr * v_interface
    end

    return SVector(f1,
                   f2 + p_interface,
                   f3, zero(eltype(u_ll)))
end

"""
	flux_tec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature1D)

Total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_tec(u_ll, u_rr, orientation::Integer,
                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    # Unpack left and right state
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
    _, _, rho_theta_ll = u_ll
    _, _, rho_theta_rr = u_rr

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = gammamean * v1_avg
    return SVector(f1, f2, f3, zero(eltype(u_ll)))
end

"""
	flux_ec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperatureWithGravity1D)

Entropy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_ec(u_ll, u_rr, orientation::Integer,
                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    # Unpack left and right state
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
    _, _, rho_theta_ll = u_ll
    _, _, rho_theta_rr = u_rr

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr) * f1
    return SVector(f1, f2, f3, zero(eltype(u_ll)))
end

"""
	flux_etec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperatureWithGravity1D)

Entropy and total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_etec(u_ll, u_rr, orientation::Integer,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    # Unpack left and right state
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    _, _, rho_theta_ll = u_ll
    _, _, rho_theta_rr = u_rr

    # Compute the necessary mean values
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f3 = gammamean * v1_avg
    f1 = f3 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    f2 = f1 * v1_avg + p_avg
    return SVector(f1, f2, f3, zero(eltype(u_ll)))
end

@inline function prim2cons(prim,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho, v1, p, phi = prim
    rho_v1 = rho * v1
    rho_theta = equations.p_0 / equations.R *
                exp(1 / equations.gamma * log(p / equations.p_0))
    return SVector(rho, rho_v1, rho_theta, phi)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho, rho_v1, rho_theta, phi = u
    v1 = rho_v1 / rho
    p = equations.K * exp(equations.gamma * log(rho_theta))
    return SVector(rho, v1, p, phi)
end

@inline function cons2cons(u,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    return u
end

@inline function cons2entropy(u,
                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho, rho_v1, rho_theta = u

    w1 = log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma
    w3 = rho / rho_theta * equations.gamma

    return SVector(w1, zero(eltype(u)), w3, zero(eltype(u)))
end

@inline function energy_total(cons,
                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    # Mathematical entropy
    p = equations.p_0 * (equations.R * cons[3] / equations.p_0)^equations.gamma

    U = p / (equations.gamma - 1) + 1 / 2 * (cons[2]^2) / (cons[1]) + cons[1] * cons[4]

    return U
end

# Default entropy is the mathematical entropy
@inline function entropy(cons,
                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    entropy_phys(cons, equations)
end

@inline function entropy_phys(cons,
                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    p = equations.K * (cons[3])^equations.gamma
    s = log(p) - equations.gamma * log(cons[1])
    S = -s * cons[1] / (equations.gamma - 1)
    return S
end

@inline function Trixi.energy_kinetic(cons,
                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    return 0.5f0 * (cons[2]^2) / (cons[1])
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho, v1, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return (abs(v1) + c,)
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate primitive variables and speed of sound
    v_mag_ll = abs(v1_ll)
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    v_mag_rr = abs(v1_rr)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    λ_max = max(v_mag_ll, v_mag_rr) + max(c_ll, c_rr)
end

@inline function pressure(cons,
                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity1D)
    p = equations.K * exp(equations.gamma * log(cons[3]))
    return p
end
end # @muladd
