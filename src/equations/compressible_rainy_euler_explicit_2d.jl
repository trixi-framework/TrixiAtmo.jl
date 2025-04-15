using Trixi
using NLsolve: nlsolve
import  Trixi: varnames,
               cons2prim, cons2entropy,
               flux, flux_chandrashekar,
               max_abs_speeds, max_abs_speed_naive,
               boundary_condition_slip_wall, entropy



###  Implementation similar to:
# Sabine Doppler, Philip L. Lederer, Joachim Schöberl, Henry von Wahl,
# A discontinuous Galerkin approach for atmospheric flows with implicit condensation,
# Journal of Computational Physics,
# Volume 499,
# 2024,
# 112713,
# ISSN 0021-9991



@muladd begin

###  equation, parameters and constants  ###

struct CompressibleRainyEulerExplicitEquations2D{RealT <: Real} <: AbstractCompressibleRainyEulerEquations{2, 7}
    # Specific heat capacities:
    c_liquid_water             ::RealT
    c_dry_air_const_pressure   ::RealT
    c_dry_air_const_volume     ::RealT
    c_vapour_const_pressure    ::RealT
    c_vapour_const_volume      ::RealT

    # Gas constants:
    R_dry_air                  ::RealT
    R_vapour                   ::RealT
    eps                        ::RealT

    # Reference values:
    ref_saturation_pressure    ::RealT
    ref_temperature            ::RealT
    ref_latent_heat_vap_temp   ::RealT
    ref_pressure               ::RealT

    # Other:
    gravity                    ::RealT
    rain_water_distr           ::RealT
    v_mean_rain                ::RealT
end


function CompressibleRainyEulerExplicitEquations2D(; RealT = Float64)
    # Specific heat capacities:
    c_liquid_water           = 4186.0
    c_dry_air_const_pressure = 1004.0
    c_dry_air_const_volume   =  717.0
    c_vapour_const_pressure  = 1885.0
    c_vapour_const_volume    = 1424.0

    # Gas constants:
    R_dry_air = c_dry_air_const_pressure - c_dry_air_const_volume
    R_vapour  = c_vapour_const_pressure  - c_vapour_const_volume
    eps       = R_dry_air                / R_vapour

    # Reference values:
    ref_saturation_pressure  = 610.7    # This needs to be adjusted if ref_temperature is changed!
    ref_temperature          = 273.15
    ref_latent_heat_vap_temp = 2.5e6#3147620.0
    ref_pressure             = 1e5

    # Other:
    gravity          = 9.81
    rain_water_distr = 8e6
    v_mean_rain      = 130.0

    return CompressibleRainyEulerExplicitEquations2D{RealT}(c_liquid_water, c_dry_air_const_pressure, c_dry_air_const_volume, 
                          c_vapour_const_pressure, c_vapour_const_volume, R_dry_air, 
                          R_vapour, eps, ref_saturation_pressure, ref_temperature,
                          ref_latent_heat_vap_temp, ref_pressure, gravity, 
                          rain_water_distr, v_mean_rain)
end



###  conversion  ###

@inline function cons2prim(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # energy density
    energy = energy_density(u, equations)

    return SVector(rho_dry, rho_vapour, rho_cloud, rho_rain, v1, v2, energy)
end


# converts consverved to entropy variables
@inline function cons2entropy(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    c_l   = equations.c_liquid_water
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # nonlinear system
    temperature = get_temperature(u, equations)
    ln_temperature  = log(temperature)
    inv_temperature = inv(temperature)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)
    v_squared_temp = 0.5 * (v1^2 + v2^2) * inv_temperature

    # check for zero density
    rho_vapour_log = 0.0

    if (rho_vapour > 0.0)
        rho_vapour_log = log(rho_vapour)
    end

    omega_dry        = c_vd * ln_temperature - R_d * log(rho_dry)   + v_squared_temp - c_vd - R_d
    omega_vapour     = c_vv * ln_temperature - R_v * rho_vapour_log + v_squared_temp - c_vv - R_v - L_ref * inv_temperature
    omega_liquid     = c_l  * ln_temperature                         + v_squared_temp - c_l
    omega_momentum_1 = -v1 * inv_temperature
    omega_momentum_2 = -v2 * inv_temperature
    omega_energy     = inv_temperature

    return SVector(omega_dry, omega_vapour, omega_liquid, omega_liquid,
                   omega_momentum_1, omega_momentum_2, omega_energy)
end


# adapted from compressible_moist_euler_2d.jl
@inline function cons2eq_pot_temp(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l      = equations.c_liquid_water
    c_pd     = equations.c_dry_air_const_pressure
    c_pv     = equations.c_vapour_const_pressure
    R_d      = equations.R_dry_air
    R_v      = equations.R_vapour
    ref_p    = equations.ref_pressure
    ref_temp = equations.ref_temperature
    ref_L    = equations.ref_latent_heat_vap_temp
    
    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # nonlinear system
    temperature = get_temperature(u, equations)

    # pressure
    p = pressure(u, equations)

    p_v  = rho_vapour * R_v * temperature
    p_d  = p - p_v
    p_vs = saturation_vapour_pressure(temperature, equations)
    H    = p_v / p_vs
    r_v  = rho_vapour / rho_dry
    r_c  = rho_cloud  / rho_dry
    r_r  = rho_rain   / rho_dry
    L_v  = ref_L + (c_pv - c_l) * temperature
    c_p  = c_pd + (r_v + r_c + r_r) * c_l

    # equivalent potential temperature
    eq_pot = (temperature * (ref_p / p_d)^(R_d / c_p) * H^(-r_v * R_v / c_p) *
               exp(L_v * r_v * inv(c_p * temperature)))

    return SVector(rho, rho_vapour, rho_cloud, rho_rain, v1, v2, eq_pot, rho)
end


# for convenience TODO rename
@inline function cons2speeds(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # get speed of sound 
    v_sound = speed_of_sound(u, equations)[1]

    # get terminal velocity rain
    v_r = terminal_velocity_rain(rho_vapour + rho_cloud, rho_rain, equations)

    return SVector(v1, v2, v_sound, v_r)
end



###  varnames  ###

varnames(::typeof(cons2cons), ::CompressibleRainyEulerExplicitEquations2D) = ("rho_dry", "rho_vapour", "rho_cloud", "rho_rain",
                                                                        "rho_v1", "rho_v2",
                                                                        "energy_density")


varnames(::typeof(cons2prim), ::CompressibleRainyEulerExplicitEquations2D) = ("rho_dry", "rho_vapour", "rho_cloud", "rho_rain",
                                                                        "v1", "v2",
                                                                        "energy_density")

varnames(::typeof(cons2eq_pot_temp), ::CompressibleRainyEulerExplicitEquations2D) = ("rho_dry", "rho_vapour",
                                                                             "rho_cloud", "rho_rain", 
                                                                             "v1", "v2", "eq_pot_temp",
                                                                             "rho")



###  physics variables  ###

@inline function densities(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # densities
    rho_dry    = u[1]
    rho_vapour = u[2]
    rho_cloud  = u[3]
    rho_rain   = u[4]
    rho        = rho_dry + rho_vapour + rho_cloud + rho_rain
    rho_inv    = inv(rho)

    return SVector(rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv)
end

@inline function rain_density(u, equations::CompressibleRainyEulerExplicitEquations2D) 
    return u[4]
end

@inline function velocities(u, rho_inv, equations::CompressibleRainyEulerExplicitEquations2D)
    return SVector(u[5] * rho_inv, u[6] * rho_inv)
end


@inline function energy_density(u, equations::CompressibleRainyEulerExplicitEquations2D)
    return u[7]
end


@inline function pressure(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    R_d = equations.R_dry_air
    R_v = equations.R_vapour

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    temperature = get_temperature(u, equations)
  
    p = (R_d * rho_dry + R_v * rho_vapour) * temperature

    return p
end


@inline function speed_of_sound(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l      = equations.c_liquid_water
    c_vd     = equations.c_dry_air_const_volume
    c_vv     = equations.c_vapour_const_volume
    R_d      = equations.R_dry_air
    R_v      = equations.R_vapour

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    if ( rho_vapour < 0.0 )
        error("rho vapour less than zero")
    end
    if ( rho_cloud < 0.0 )
        error("rho cloud less than zero")
    end

    # formula
    p       = pressure(u, equations)
    q_v     = rho_vapour / rho_dry
    q_l     = (rho_cloud + rho_rain) / rho_dry
    gamma_m = 1 + (R_d + R_v * q_v) / (c_vd + c_vv * q_v + c_l * q_l)

    if (rho_inv < 0.0)
        error("rho less than zero")
    elseif (p < 0.0)
        error("pressure less than zero")
    end

    v_sound = sqrt(gamma_m * p * rho_inv)
    
    return SVector(v_sound, gamma_m)
end


@inline function get_temperature(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l   = equations.c_liquid_water
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    L_ref = equations.ref_latent_heat_vap_temp

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)

    # energy density
    energy = energy_density(u, equations::CompressibleRainyEulerExplicitEquations2D)

    return (energy - L_ref * rho_vapour - 0.5 * rho * (v1^2 + v2^2)) / (c_vd * rho_dry + c_vv * rho_vapour + c_l * (rho_cloud + rho_rain))
end


@inline function terminal_velocity_rain(rho_moist, rho_rain, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    N_0 = equations.rain_water_distr
    v_0 = equations.v_mean_rain

    # formula ( \Gamma(4.5) / 6 ~= 1.9386213994279082 )
    if ( rho_rain > 0.0)
        v_terminal_rain = v_0 * 1.9386213994279082 * (rho_rain / (pi * (rho_moist + rho_rain) * N_0))^(0.125)
    else
        v_terminal_rain = 0.0
    end

    return v_terminal_rain
end


@inline function saturation_vapour_pressure(temperature, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l      = equations.c_liquid_water
    c_pv     = equations.c_vapour_const_pressure
    R_v      = equations.R_vapour
    ref_s_p  = equations.ref_saturation_pressure
    ref_temp = equations.ref_temperature
    ref_L    = equations.ref_latent_heat_vap_temp

    # testing 
    if (temperature < 0.0)
        display(temperature)
        error("temp less than zero")
    end

    # Clausius Clapeyron formula
    p_vapour_saturation  = ref_s_p * (temperature / ref_temp)^((c_pv - c_l) / R_v)
    p_vapour_saturation *= exp(((ref_L - (c_pv - c_l) * ref_temp) / R_v) * (1 / ref_temp - 1 / temperature))

    return p_vapour_saturation
end


@inline function saturation_vapour_pressure_derivative(temperature, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l      = equations.c_liquid_water
    c_pv     = equations.c_vapour_const_pressure
    R_v      = equations.R_vapour
    ref_s_p  = equations.ref_saturation_pressure
    ref_temp = equations.ref_temperature
    ref_L    = equations.ref_latent_heat_vap_temp

    # testing 
    if (temperature < 0.0)
        display(temperature)
        error("temp less than zero")
    end

    const_1 = (c_pv - c_l) / R_v
    const_2 = (ref_L - (c_pv - c_l) * ref_temp) / R_v

    p_vapour_saturation_derivative  = ref_s_p / (ref_temp^const_1)
    p_vapour_saturation_derivative *= (const_1 * temperature^(const_1 - 1) + const_2 * temperature^(const_1 - 2))
    p_vapour_saturation_derivative *= exp(const_2 * (1 / ref_temp - 1 / temperature))

    return p_vapour_saturation_derivative
end


# adapted from compressible_moist_euler_2d.jl
@inline function moist_air_phase_change(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    R_v = equations.R_vapour

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # temperature
    temperature = get_temperature(u, equations)

    # saturation vapor pressure
    p_vs = saturation_vapour_pressure(temperature, equations)

    # saturation density of vapor
    rho_star_qv = p_vs / (R_v * temperature)

    # Fisher-Burgmeister-Function
    a = rho_star_qv - rho_vapour
    b = rho_cloud

    # saturation control factor
    # < 1: stronger saturation effect
    # > 1: weaker saturation effect
    C = 1.0

    return (a + b - sqrt(a^2 + b^2)) * C
end



###  pde discretization  ###

@inline function flux(u, orientation::Integer, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l = equations.c_liquid_water

    #densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # recover rho_vapour, rho_cloud, temperature from nonlinear system
    temperature = get_temperature(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)
    v_r    = terminal_velocity_rain(rho_vapour + rho_cloud, rho_rain, equations)

    # pressure
    p = pressure(u, equations)

    # energy density
    energy = energy_density(u, equations)

    # flux for orientation cases 
    if (orientation == 1)
        # "mass"
        f1 = rho_dry    * v1
        f2 = rho_vapour * v1
        f3 = rho_cloud  * v1
        f4 = rho_rain   * v1

        # "momentum"
        f5 = rho * v1 * v1 + p
        f6 = rho * v1 * v2

        # "energy"
        f7 = (energy + p) * v1

    else
        # "mass"
        f1 = rho_dry    * v2
        f2 = rho_vapour * v2
        f3 = rho_cloud  * v2
        f4 = rho_rain   * (v2 - v_r)

        # "momentum"
        f5 = rho  * v1 * v2      - rho_rain * v_r * v1
        f6 = rho  * v2 * v2  + p - rho_rain * v_r * v2

        # "energy"
        f7 = (energy + p) * v2 - (c_l * temperature + 0.5 * (v1^2 + v2^2)) * rho_rain * v_r
    end

    return SVector(f1, f2, f3, f4, f5, f6, f7)
end


@inline function flux(u, normal_direction::AbstractVector, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l = equations.c_liquid_water

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # recover rho_vapour, rho_cloud, temperature from nonlinear system
    temperature = get_temperature(u, equations)

    # velocity
    v1, v2 = velocities(u, rho_inv, equations)
    v_r    = terminal_velocity_rain(rho_vapour + rho_cloud, rho_rain, equations)

    # normal velocities
    v_normal   = v1  * normal_direction[1] +  v2 * normal_direction[2]
    v_r_normal =                             v_r * normal_direction[2]

    # pressure
    p = pressure(u, equations)

    # energy density
    energy = energy_density(u, equations)

    # flux
    # "mass"
    f1 = rho_dry    *  v_normal
    f2 = rho_vapour *  v_normal
    f3 = rho_cloud  *  v_normal
    f4 = rho_rain   * (v_normal - v_r_normal)

    # "momentum"
    f5 = rho * v_normal * v1 + p * normal_direction[1] - rho_rain * v_r_normal * v1
    f6 = rho * v_normal * v2 + p * normal_direction[2] - rho_rain * v_r_normal * v2 

    # "energy"
    f7 = (energy + p) * v_normal - (c_l * temperature + 0.5 * (v1^2 + v2^2)) * rho_rain * v_r_normal

    return SVector(f1, f2, f3, f4, f5, f6, f7)
end


# no Coriolis term
@inline function source_terms_rainy(u, x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    R_v      = equations.R_vapour
    ref_temp = equations.ref_temperature
    g        = equations.gravity
    
    # name needed variables
    rho_v2 = u[6]

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    # temperature
    temperature = get_temperature(u, equations)

    rho_vs = saturation_vapour_pressure(temperature, equations) / (R_v * temperature)
    Q_ph   = moist_air_phase_change(u, equations)

    # source terms phase change
    S_evaporation     = (3.86e-3 - 9.41e-5 * (temperature - ref_temp)) * (1 + 9.1 * rho_rain^(0.1875))
    S_evaporation    *= (rho_vs - rho_vapour) * rho_rain^(0.5)
    S_auto_conversion = 0.001 * rho_cloud
    S_accretion       = 1.72 * rho_cloud * rho_rain^(0.875)
    S_rain            = S_auto_conversion + S_accretion - S_evaporation

    return SVector(0.0, Q_ph + S_evaporation, -Q_ph - S_auto_conversion - S_accretion, S_rain, 0.0,
                   -rho * g, -rho_v2 * g)
end


# no Coriolis term
@inline function source_terms_moist(u, x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    g = equations.gravity
    
    # name needed variables
    rho_v2 = u[6]

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    Q_ph = moist_air_phase_change(u, equations)

    return SVector(0.0, Q_ph, -Q_ph, 0.0, 0.0,
                   -rho * g, -rho_v2 * g)
end


# no phase changes and no Coriolis term
@inline function source_terms_no_phase_change(u, x, t, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    g = equations.gravity
    
    # name needed variables
    rho_v2 = u[6]

    # densities
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)

    return SVector(0.0, 0.0, 0.0, 0.0, 0.0, -g * rho, -g * rho_v2) 
end


@inline function max_abs_speeds(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # name needed variables
    v1, v2, v_sound, v_r = cons2speeds(u, equations)

    return SVector((abs(v1) + v_sound), (abs(v2) + v_sound + abs(v_r)))
end


@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleRainyEulerExplicitEquations2D)
    # name needed variables
    v1_ll, v2_ll, v_sound_ll, v_r_ll = cons2speeds(u_ll, equations)
    v1_rr, v2_rr, v_sound_rr, v_r_rr = cons2speeds(u_rr, equations)

    # calculate upper bounds for left and right speed
    v_ll_max  = abs(v1_ll * normal_direction[1] +  v2_ll * normal_direction[2])
    v_ll_max += abs(                              v_r_ll * normal_direction[2])

    v_rr_max  = abs(v1_rr * normal_direction[1] +  v2_rr * normal_direction[2])
    v_rr_max += abs(                              v_r_rr * normal_direction[2])

    return max(v_ll_max, v_rr_max) + max(v_sound_ll, v_sound_rr) * norm(normal_direction)
end


@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::CompressibleRainyEulerExplicitEquations2D)
    # name needed variables
    v1_ll, v2_ll, v_sound_ll, v_r_ll = cons2speeds(u_ll, equations)
    v1_rr, v2_rr, v_sound_rr, v_r_rr = cons2speeds(u_rr, equations)

    if (orientation == 1)
        v_ll  = abs(v1_ll)
        v_rr  = abs(v1_rr)
    else
        v_ll  = abs(v2_ll)
        v_rr  = abs(v2_rr)
    end
    # experimental
    return max(v_ll, v_rr) + max(v_sound_ll, v_sound_rr) + max(abs(v_r_ll), abs(v_r_rr))
end


###  boundary conditions  ###

# adapted from compressible_moist_euler_2d.jl which was probably adapted from compressible_euler_2d.jl
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector, x, t, 
                                              surface_flux_function, equations::CompressibleRainyEulerExplicitEquations2D)

    norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    # rotate the internal solution state
    u_local = rotate_to_x(u_inner, normal, equations)
    
    # name needed variables
    rho_v1 = u_local[5]

    # densities
    rho_dry_local, rho_vapour_local, rho_cloud_local, rho_rain_local, rho_local, rho_inv_local = densities(u_local, equations)

    # velocities
    v_normal       = rho_v1 * rho_inv_local
    v_sound, gamma = speed_of_sound(u_local, equations)

    # pressure
    p_local = pressure(u_local, equations)
    
    # Get the solution of the pressure Riemann problem
    # See Section 6.3.3 of
    # Eleuterio F. Toro (2009)
    # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
    # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
    if v_normal <= 0.0
        p_star = p_local *
                 (1.0 + 0.5 * (gamma - 1) * v_normal / v_sound)^(2.0 * gamma *
                                                                     inv(gamma - 1))
    else # v_normal > 0.0
        A = 2.0 / ((gamma + 1) * rho_local)
        B = p_local * (gamma - 1) / (gamma + 1)
        p_star = p_local +
                 0.5 * v_normal / A *
                 (v_normal + sqrt(v_normal^2 + 4.0 * A * (p_local + B)))
    end

    # For the slip wall we directly set the flux as the normal velocity is zero
    return SVector(0.0, 0.0, 0.0, 0.0,
                   p_star * normal[1] * norm_,
                   p_star * normal[2] * norm_,
                   0.0)
end


# same as in compressible_euler_2d.jl
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector, direction, x, t,
                                              surface_flux_function, equations::CompressibleRainyEulerExplicitEquations2D)
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


# same as in compressible_euler_2d.jl
@inline function rotate_to_x(u, normal_vector, equations::CompressibleRainyEulerExplicitEquations2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    return SVector(u[1], u[2], u[3], u[4],
                   c * u[5] + s * u[6],
                  -s * u[5] + c * u[6],
                   u[7])
end


# should be used together with TreeMesh (adapted from compressible_euler_2d.jl)
@inline function boundary_condition_slip_wall(u_inner, orientation, direction, x, t, surface_flux_function,
                                              equations::CompressibleRainyEulerExplicitEquations2D)
    # get the appropriate normal vector from the orientation
    RealT = eltype(u_inner)
    if orientation == 1
        normal_direction = SVector(one(RealT), zero(RealT))
    else # orientation == 2
        normal_direction = SVector(zero(RealT), one(RealT))
    end

    # compute and return the flux using `boundary_condition_slip_wall` routine above
    return boundary_condition_slip_wall(u_inner, normal_direction, direction,
                                        x, t, surface_flux_function, equations)
end


#= for parabolic terms (LaplaceDiffusion2D)
@inline function boundary_condition_laplace(flux_inner, u_inner, normal::AbstractVector, x, t, operator_type::Trixi.Gradient,
                                            equations_parabolic::LaplaceDiffusion2D)
    return u_inner
end

@inline function boundary_condition_laplace(flux_inner, u_inner, normal::AbstractVector, x, t, operator_type::Trixi.Divergence,
                                            equations_parabolic::LaplaceDiffusion2D)
    return flux_inner
end=#


# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.
# adapted from compressible_moist_euler_2d.jl, does NOT work with rain!
@inline function flux_LMARS(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    a = 360.0

    # densities
    rho_dry_ll, rho_vapour_ll, rho_cloud_ll, rho_rain_ll, rho_ll, rho_inv_ll = densities(u_ll, equations)
    rho_dry_rr, rho_vapour_rr, rho_cloud_rr, rho_rain_rr, rho_rr, rho_inv_rr = densities(u_rr, equations)

    # pressure
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    # velocities
    v1_ll, v2_ll = velocities(u_ll, rho_inv_ll, equations)
    v1_rr, v2_rr = velocities(u_rr, rho_inv_rr, equations)
    
    v_ll  = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr  = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    norm_ = norm(normal_direction)

    # diffusion parameter 0.0 < beta <= 1.0
    beta = 1.0
    
    # interface flux components
    rho = 0.5 * (rho_ll + rho_rr)
    p_interface = 0.5 * (p_ll + p_rr) - beta * 0.5 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5 * (v_ll + v_rr) - beta * 1 / (2 * a * rho) * (p_rr - p_ll) * norm_
    
    if (v_interface > 0)
        f1, f2, f3, _, f5, f6, f7 = u_ll * v_interface
        f7 += p_ll * v_interface
    else
        f1, f2, f3, _, f5, f6, f7 = u_rr * v_interface
        f7 += p_rr * v_interface
    end
    
    return SVector(f1, f2, f3, 0.0,
                   f5 + p_interface * normal_direction[1],
                   f6 + p_interface * normal_direction[2],
                   f7)
end


# Adjusted EC flux in a normal direction with R_q=0. This is based on
# A. Gouasmi, K. Duraisamy, S. M. Murman, Formulation of Entropy-Stable schemes for the
# multicomponent compressible Euler equations, 4 Feb 2020, doi:10.1016/j.cma.2020.112912,
# https://arxiv.org/abs/1904.00972 [math.NA].
# Adapted from compressible_moist_euler_2d.jl, careful, does NOT work for rain!
@inline function flux_chandrashekar(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l   = equations.c_liquid_water
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    ref_L = equations.ref_latent_heat_vap_temp
    R_q   = 0.0

    # densities and temperatures
    rho_dry_ll, rho_vapour_ll, rho_cloud_ll, rho_rain_ll, rho_ll, rho_inv_ll = densities(u_ll, equations)
    rho_dry_rr, rho_vapour_rr, rho_cloud_rr, rho_rain_rr, rho_rr, rho_inv_rr = densities(u_rr, equations)
    temperature_ll                                                           = get_temperature(u_ll, equations)
    temperature_rr                                                           = get_temperature(u_rr, equations)

    # velocities
    v1_ll, v2_ll = velocities(u_ll, rho_inv_ll, equations)
    v1_rr, v2_rr = velocities(u_rr, rho_inv_rr, equations)
    vr_ll        = terminal_velocity_rain(rho_vapour_ll + rho_cloud_ll, rho_rain_ll, equations)
    vr_rr        = terminal_velocity_rain(rho_vapour_rr + rho_cloud_rr, rho_rain_rr, equations)
    
    # mean values
    rho_dry_mean          = 0.0
    rho_vapour_mean       = 0.0
    rho_cloud_mean        = 0.0
    rho_rain_mean         = 0.0
    inv_temperature_mean  = 0.0

    if (!(rho_dry_ll == 0.0) && !(rho_dry_rr == 0.0))
        rho_dry_mean = ln_mean(rho_dry_ll, rho_dry_rr)
    end

    if (!(rho_vapour_ll == 0.0) && !(rho_vapour_rr == 0.0))
        rho_vapour_mean = ln_mean(rho_vapour_ll, rho_vapour_rr)
    end

    if (!(rho_cloud_ll == 0.0) && !(rho_cloud_rr == 0.0))
        rho_cloud_mean = ln_mean(rho_cloud_ll, rho_cloud_rr)
    end

    if (!(rho_rain_ll == 0.0) && !(rho_rain_rr == 0.0))
        rho_rain_mean = ln_mean(rho_rain_ll, rho_rain_rr)
    end

    if (!(inv(temperature_ll) == 0.0) && !(inv(temperature_rr) == 0.0))
        inv_temperature_mean = inv_ln_mean(inv(temperature_ll), inv(temperature_rr))
    end
    
    v1_avg              = 0.5 * (v1_ll                 + v1_rr)
    v2_avg              = 0.5 * (v2_ll                 + v2_rr)
    v1_square_avg       = 0.5 * (v1_ll^2               + v1_rr^2)
    v2_square_avg       = 0.5 * (v2_ll^2               + v2_rr^2)
    rho_dry_avg         = 0.5 * (rho_dry_ll            + rho_dry_rr)
    rho_vapour_avg      = 0.5 * (rho_vapour_ll         + rho_vapour_rr)
    rho_cloud_avg       = 0.5 * (rho_cloud_ll          + rho_cloud_rr)
    rho_rain_avg        = 0.5 * (rho_rain_ll           + rho_rain_rr)
    inv_temperature_avg = 0.5 * (inv(temperature_ll)   + inv(temperature_rr))
    v_dot_n_avg         = normal_direction[1] * v1_avg + normal_direction[2] * v2_avg
    
    p_int  = inv(inv_temperature_avg) * (R_d * rho_dry_avg + R_v * rho_vapour_avg + R_q * (rho_cloud_avg + rho_rain_avg))
    K_avg  = 0.5 * (v1_square_avg + v2_square_avg)
    vr_normal_avg = 0.5 * normal_direction[2] * (vr_ll + vr_rr)
    rho_rain_avg        = 0.5 * (rho_rain_ll   + rho_rain_rr)

    # assemble the flux
    f_dry    = rho_dry_mean    *  v_dot_n_avg
    f_vapour = rho_vapour_mean *  v_dot_n_avg
    f_cloud  = rho_cloud_mean  *  v_dot_n_avg
    f_rain   = rho_rain_avg    * (v_dot_n_avg - vr_normal_avg)
    f_rhov1  = (f_dry + f_vapour + f_cloud + f_rain) * v1_avg + normal_direction[1] * p_int
    f_rhov2  = (f_dry + f_vapour + f_cloud + f_rain) * v2_avg + normal_direction[2] * p_int
    f_energy = ((        c_vd * inv_temperature_mean - K_avg)  *  f_dry    +
                (ref_L + c_vv * inv_temperature_mean - K_avg)  *  f_vapour +
                (         c_l * inv_temperature_mean - K_avg)  * (f_cloud  + f_rain)  + v1_avg * f_rhov1 + v2_avg * f_rhov2)

    return SVector(f_dry, f_vapour, f_cloud, f_rain, f_rhov1, f_rhov2, f_energy)
end


@inline function flux_chandrashekar(u_ll, u_rr, orientation::Int, equations::CompressibleRainyEulerExplicitEquations2D)
    if (orientation == 1)
        return flux_chandrashekar(u_ll, u_rr, SVector(1, 0), equations)
    else
        return flux_chandrashekar(u_ll, u_rr, SVector(0, 1), equations)
    end
end



@inline function flux_ec_rain(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l   = equations.c_liquid_water
    c_vd  = equations.c_dry_air_const_volume
    c_vv  = equations.c_vapour_const_volume
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    L_ref = equations.ref_latent_heat_vap_temp

    # densities and temperatures
    rho_dry_ll, rho_vapour_ll, rho_cloud_ll, rho_rain_ll, rho_ll, rho_inv_ll = densities(u_ll, equations)
    rho_dry_rr, rho_vapour_rr, rho_cloud_rr, rho_rain_rr, rho_rr, rho_inv_rr = densities(u_rr, equations)
    temperature_ll                                                           = get_temperature(u_ll, equations)
    temperature_rr                                                           = get_temperature(u_rr, equations)
    inv_temperature_ll                                                       = inv(temperature_ll)
    inv_temperature_rr                                                       = inv(temperature_rr)

    # velocities
    v1_ll, v2_ll = velocities(u_ll, rho_inv_ll, equations)
    v1_rr, v2_rr = velocities(u_rr, rho_inv_rr, equations)
    vr_ll        = terminal_velocity_rain(rho_vapour_ll + rho_cloud_ll, rho_rain_ll, equations)
    vr_rr        = terminal_velocity_rain(rho_vapour_rr + rho_cloud_rr, rho_rain_rr, equations)

    # velocity averages
    v1_avg        = 0.5 * (v1_ll                 + v1_rr)
    v2_avg        = 0.5 * (v2_ll                 + v2_rr)
    v1_square_avg = 0.5 * (v1_ll^2               + v1_rr^2)
    v2_square_avg = 0.5 * (v2_ll^2               + v2_rr^2)
    K_avg         = 0.5 * (v1_square_avg         + v2_square_avg)
    v_dot_n_avg   = normal_direction[1] * v1_avg + normal_direction[2] *  v2_avg
    vr_normal_avg =                          0.5 * normal_direction[2] * (vr_ll + vr_rr)

    # density averages
    rho_dry_avg         = 0.5 * (rho_dry_ll    + rho_dry_rr)
    rho_vapour_avg      = 0.5 * (rho_vapour_ll + rho_vapour_rr)
    rho_cloud_avg       = 0.5 * (rho_cloud_ll  + rho_cloud_rr)
    rho_rain_avg        = 0.5 * (rho_rain_ll   + rho_rain_rr)

    # density log means
    rho_dry_log          = 0.0
    rho_vapour_log       = 0.0
    rho_cloud_log        = 0.0
    rho_rain_log         = 0.0

    if (!(rho_dry_ll == 0.0) && !(rho_dry_rr == 0.0))
        rho_dry_log = ln_mean(rho_dry_ll, rho_dry_rr)
    end

    if (!(rho_vapour_ll == 0.0) && !(rho_vapour_rr == 0.0))
        rho_vapour_log = ln_mean(rho_vapour_ll, rho_vapour_rr)
    end

    if (!(rho_cloud_ll == 0.0) && !(rho_cloud_rr == 0.0))
        rho_cloud_log = ln_mean(rho_cloud_ll, rho_cloud_rr)
    end

    if (!(rho_rain_ll == 0.0) && !(rho_rain_rr == 0.0))
        rho_rain_log = ln_mean(rho_rain_ll, rho_rain_rr)
    end

    # other averages
    inv_temperature_avg = 0.5 * (inv_temperature_ll + inv_temperature_rr)
    inv_temperature_log = inv_ln_mean(inv_temperature_ll, inv_temperature_rr)
    p_int  = inv(inv_temperature_avg) * (R_d * rho_dry_avg + R_v * rho_vapour_avg)
    
    # density flux
    f_vapour = rho_vapour_log * v_dot_n_avg
    f_cloud  = rho_cloud_avg  * v_dot_n_avg
    f_dry    = rho_dry_log  *  v_dot_n_avg
    f_rain   = rho_rain_avg * (v_dot_n_avg - vr_normal_avg)
    f_moist  = f_vapour + f_cloud
    f_rho    = f_dry    + f_moist + f_rain

    # momentum flux
    f_rhov1  = f_rho * v1_avg + p_int * normal_direction[1]
    f_rhov2  = f_rho * v2_avg + p_int * normal_direction[2]

    # energy flux
    f_energy = (c_vd * inv_temperature_log - K_avg) * f_dry + (c_vv * inv_temperature_log - K_avg + L_ref) * f_vapour + 
               (c_l  * inv_temperature_log - K_avg) * (f_cloud + f_rain) + (v1_avg * f_rhov1 + v2_avg * f_rhov2) 

    return SVector(f_dry, f_vapour, f_cloud, f_rain, f_rhov1, f_rhov2, f_energy)
end

@inline function flux_ec_rain(u_ll, u_rr, orientation::Int, equations::CompressibleRainyEulerExplicitEquations2D)
    if (orientation == 1)
        return flux_ec_rain(u_ll, u_rr, SVector(1, 0), equations)
    else
        return flux_ec_rain(u_ll, u_rr, SVector(0, 1), equations)
    end
end



# adapted from ShallowWaterEquations2D (Recommended with rain!)
@inline function boundary_condition_simple_slip_wall(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function,
                                                     equations::CompressibleRainyEulerExplicitEquations2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
    u_normal = normal[1] * u_inner[5] + normal[2] * u_inner[6]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], u_inner[4],
    u_inner[5] - 2 * u_normal * normal[1],
    u_inner[6] - 2 * u_normal * normal[2],
    u_inner[7])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

    return flux
end


# adapted from ShallowWaterEquations2D (Recommended with rain!)
@inline function boundary_condition_simple_slip_wall(u_inner, orientation, direction, x, t, surface_flux_function,
                                                     equations::CompressibleRainyEulerExplicitEquations2D)
    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], u_inner[4], -u_inner[5],  u_inner[6], u_inner[7])
    else # orientation == 2
        u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3],  u_inner[4], u_inner[5], -u_inner[6], u_inner[7])
    end

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end


@inline function entropy(u, equations::CompressibleRainyEulerExplicitEquations2D)
    # constants
    c_l  = equations.c_liquid_water
    c_vd = equations.c_dry_air_const_volume
    c_vv = equations.c_vapour_const_volume
    R_d  = equations.R_dry_air
    R_v  = equations.R_vapour

    # variables
    rho_dry, rho_vapour, rho_cloud, rho_rain, rho, rho_inv = densities(u, equations)
    temperature = get_temperature(u, equations)

    # s_k
    s_d = c_vd * log(temperature) - R_d * log(rho_dry)
    s_v = c_vv * log(temperature) - R_v * log(rho_vapour)
    s_l = c_l  * log(temperature)
    
    return rho_dry * s_d + rho_vapour * s_v + (rho_cloud + rho_rain) * s_l
end

end  # muladd end