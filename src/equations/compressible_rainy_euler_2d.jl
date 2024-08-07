using Trixi
import Trixi: varnames,
              cons2prim, cons2entropy,
              flux,
              max_abs_speeds, max_abs_speed_naive,
              boundary_condition_slip_wall



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

struct CompressibleRainyEulerEquations2D{RealT <: Real} <: AbstractCompressibleMoistEulerEquations{2, 15}
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


function CompressibleRainyEulerEquations2D(; RealT = Float64)
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
    saturation_vapour_pressure = 610.7
    ref_temperature            = 273.15
    latent_heat_vap_ref_temp   = 2.5e6
    ref_pressure               = 1e5

    # Other:
    gravity          = 9.81
    rain_water_distr = 8e6
    v_mean_rain      = 130.0

    return CompressibleRainyEulerEquations2D{RealT}(c_liquid_water, c_dry_air_const_pressure, c_dry_air_const_volume, 
                          c_vapour_const_pressure, c_vapour_const_volume, R_dry_air, 
                          R_vapour, eps, saturation_vapour_pressure, ref_temperature,
                          latent_heat_vap_ref_temp, ref_pressure, gravity, 
                          rain_water_distr, v_mean_rain)
end



###  varnames  ###

varnames(::typeof(cons2cons), ::CompressibleRainyEulerEquations2D) = ("rho_dry_", "rho_moist_", "rho_rain_",
                                                                      "rho_v1", "rho_v2",
                                                                      "energy_",
                                                                      "rho_vapour_h", "rho_cloud_h",
                                                                      "temperature",
                                                                      "rho_dry_h_0", "rho_moist_h_0", "rho_rain_h_0",
                                                                      "energy_h_0",
                                                                      "rho_vapour_h_0",
                                                                      "temperature_0")


varnames(::typeof(cons2prim), ::CompressibleRainyEulerEquations2D) = ("rho_dry", "rho_moist", "rho_rain",
                                                                      "v1", "v2",
                                                                      "energy",
                                                                      "rho_vapour", "rho_cloud",
                                                                      "temperature",
                                                                      "rho_dry_h_0", "rho_moist_h_0", "rho_rain_h_0",
                                                                      "energy_h_0",
                                                                      "rho_vapour_h_0",
                                                                      "temperature_0")




###  conversion  ###

@inline function cons2prim(u, equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    rho_v1, rho_v2,
    energy_,
    _, _, _,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    energy_h_0,
    _, _      = u

    # densities
    rho_dry   = rho_dry_     + rho_dry_h_0
    rho_moist = rho_moist_   + rho_moist_h_0
    rho_rain  = rho_rain_    + rho_rain_h_0
    rho_inv   = 1 / (rho_dry + rho_moist + rho_rain)

    # velocity
    v1        = rho_v1       * rho_inv
    v2        = rho_v2       * rho_inv

    # energy
    energy    = energy_      + energy_h_0

    return SVector(rho_dry, rho_moist, rho_rain, v1, v2, energy,
                   u[7], u[8], u[9], u[10], u[11], u[12], u[13], u[14], u[15])
end


@inline function cons2entropy(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end


# for convenience
@inline function cons2speeds(u, equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    rho_v1, rho_v2,
    _, _, _, _,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    _, _, _   = u

    # densities
    rho_dry   = rho_dry_     + rho_dry_h_0
    rho_moist = rho_moist_   + rho_moist_h_0
    rho_rain  = rho_rain_    + rho_rain_h_0
    rho_inv   = 1 / (rho_dry + rho_moist + rho_rain)

    # velocity
    v1        = rho_v1       * rho_inv
    v2        = rho_v2       * rho_inv

    # get speed of sound 
    v_sound   = speed_of_sound(u, equations)

    # get terminal velocity rain
    v_r       = terminal_velocity_rain(rho_moist, rho_rain, equations)

    return SVector(v1, v2, v_sound, v_r)
end



###  physics variables  ###

@inline function speed_of_sound(u, equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l      = equations.c_liquid_water
    c_vd     = equations.c_dry_air_const_volume
    c_vv     = equations.c_vapour_const_volume
    R_d      = equations.R_dry_air
    R_v      = equations.R_vapour

    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    _, _,
    _,
    rho_vapour_h, rho_cloud_h,
    temperature,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    _,
    _, _       = u

    # densities
    rho_dry    = rho_dry_     + rho_dry_h_0
    rho_moist  = rho_moist_   + rho_moist_h_0
    rho_rain   = rho_rain_    + rho_rain_h_0
    rho_vapour = rho_vapour_h     #TODO double check formula because this seems weird at first
    rho_cloud  = rho_cloud_h      #TODO ^
    rho_inv    = 1 / (rho_dry + rho_moist + rho_rain)

    # formula
    p   = R_d * rho_dry * temperature + R_v * rho_vapour * temperature
    q_v = rho_vapour / rho_dry
    q_l = (rho_cloud + rho_rain) / rho_dry
    gamma_m = 1 + (R_d + R_v * q_v) / (c_vd + c_vv * q_v + c_l * q_l)
    v_sound = sqrt(gamma_m * p * rho_inv)
    
    return v_sound
end


@inline function terminal_velocity_rain(rho_moist, rho_rain, equations::CompressibleRainyEulerEquations2D)
    # constants
    N_0 = equations.rain_water_distr
    v_0 = equations.v_mean_rain

    # formula ( gamma(4.5) / 6 ~= 1.9386213994279082 )
    v_terminal_rain = v_0 * 1.9386213994279082 * (rho_rain / (pi * (rho_moist + rho_rain) * N_0))^(0.125)

    return v_terminal_rain
end


@inline function saturation_vapour_pressure(temperature, equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l      = equations.c_liquid_water
    c_pv     = equations.c_vapour_const_pressure
    R_v      = equations.R_vapour
    ref_s_p  = equations.ref_saturation_pressure
    ref_temp = equations.ref_temperature
    ref_L    = equations.ref_latent_heat_vap_temp

    # Clausius Clapeyron formula
    p_vapour_saturation  = ref_s_p * (temperature / ref_temp)^((c_pv - c_l) / R_v)
    p_vapour_saturation *= exp(((ref_L - (c_pv - c_l) * ref_temp) / R_v) * (1 / ref_temp - 1 / temperature))

    return p_vapour_saturation
end



###  pde discretization  ###

@inline function flux(u, orientation::Integer, equations::CompressibleRainyEulerEquations2D)
    # constants
    c_l   = equations.c_liquid_water
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    T_ref = equations.ref_temperature

    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    rho_v1, rho_v2,
    energy_,
    rho_vapour_h, _, 
    temperature,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    energy_h_0,
    rho_vapour_h_0, 
    temperature_0  = u

    # densities
    rho_dry    = rho_dry_     + rho_dry_h_0
    rho_moist  = rho_moist_   + rho_moist_h_0
    rho_rain   = rho_rain_    + rho_rain_h_0
    rho_vapour = rho_vapour_h                              #TODO ?
    rho        = rho_dry      + rho_moist      + rho_rain
    rho_inv    = 1 / (rho)

    # velocity
    v1         = rho_v1       * rho_inv
    v2         = rho_v2       * rho_inv
    v_r        = terminal_velocity_rain(rho_moist, rho_rain, equations)

    # pressure
    p   = R_d * rho_dry     * temperature   + R_v * rho_vapour     * temperature
    p_0 = R_d * rho_dry_h_0 * temperature_0 + R_v * rho_vapour_h_0 * temperature_0
    p_  = p   - p_0

    # energy
    energy = energy_ + energy_h_0

    # flux for orientation cases 
    if (orientation == 1)
        # "mass"
        f1 = rho_dry   * v1
        f2 = rho_moist * v1
        f3 = rho_rain  * v1

        # "impulse"
        f4 = rho       * v1 * v1 + p_
        f5 = rho       * v1 * v2

        # "energy"
        f6 = (energy + p) * v1

    else
        # "mass"
        f1 = rho_dry   *  v2
        f2 = rho_moist *  v2
        f3 = rho_rain  * (v2 - v_r)

        # "impulse"
        f4 = rho       * v1 * v2 - rho_rain * v_r * v1
        f5 = rho       * v2 * v2 - rho_rain * v_r * v2 + p_

        # "energy"
        f6 = ((energy + p) * v2 - (c_l * (temperature - T_ref) + 0.5 * (v1^2 + v2^2))) * rho_rain * v_r
    end

    return SVector(f1, f2, f3, f4, f5, f6, 
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end


@inline function flux(u, normal_direction::AbstractVector, equations::CompressibleRainyEulerEquations2D)
    #TODO Double-check for mistakes in "impulse" and "energy"!
    # constants
    c_l   = equations.c_liquid_water
    R_d   = equations.R_dry_air
    R_v   = equations.R_vapour
    T_ref = equations.ref_temperature

    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    rho_v1, rho_v2,
    energy_,
    rho_vapour_h, _, 
    temperature,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    energy_h_0,
    rho_vapour_h_0, 
    temperature_0  = u

    # densities
    rho_dry    = rho_dry_     + rho_dry_h_0
    rho_moist  = rho_moist_   + rho_moist_h_0
    rho_rain   = rho_rain_    + rho_rain_h_0
    rho_vapour = rho_vapour_h                              #TODO ?
    rho        = rho_dry      + rho_moist      + rho_rain
    rho_inv    = 1 / (rho)

    # velocity
    v1         = rho_v1       * rho_inv
    v2         = rho_v2       * rho_inv
    v_r        = terminal_velocity_rain(rho_moist, rho_rain, equations)

    # normal velocities
    v_normal   = v1  * normal_direction[1] +  v2 * normal_direction[2]
    v_r_normal =                             v_r * normal_direction[2]    #TODO correct?

    # pressure
    p   = R_d * rho_dry     * temperature   + R_v * rho_vapour     * temperature
    p_0 = R_d * rho_dry_h_0 * temperature_0 + R_v * rho_vapour_h_0 * temperature_0
    p_  = p   - p_0

    # energy
    energy = energy_ + energy_h_0

    # flux
    # "mass"
    f1 = rho_dry   *  v_normal
    f2 = rho_moist *  v_normal
    f3 = rho_rain  * (v_normal - v_r_normal)

    # "impulse"
    f4 = rho       * v_normal * v1 - rho_rain * v_r_normal * v1 + p_ * normal_direction[1]
    f5 = rho       * v_normal * v2 - rho_rain * v_r_normal * v2 + p_ * normal_direction[2]

    # "energy"
    f6 = ((energy + p) * v_normal - (c_l * (temperature - T_ref) + 0.5 * (v_normal^2))) * rho_rain * v_r_normal

    return SVector(f1, f2, f3, f4, f5, f6, 
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end


@inline function source_terms_rainy(u, x, t, equations::CompressibleRainyEulerEquations2D)
    #TODO implement conversion source terms
    # constants
    g = equations.gravity
    
    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    _, rho_v2,
    _,
    rho_vapour_h, rho_cloud_h,
    _,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    _,
    _, _      = u

    # densities
    rho_dry   = rho_dry_    + rho_dry_h_0
    rho_moist = rho_moist_  + rho_moist_h_0
    rho_rain  = rho_rain_   + rho_rain_h_0
    
    rho       = rho_dry     + rho_moist     + rho_rain
    rho_h_0   = rho_dry_h_0 + rho_moist_h_0 + rho_rain_h_0
    rho_      = rho         - rho_h_0

    #TODO non-linear solve for temperature, rho_vapour_h, rho_cloud_h

    return SVector(0.0, 0.0, 0.0, 0.0,
                   rho_ * g, 0.0, rho_v2 * g,     #TODO check signs
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end


# no phase changes and no Coriolis term

@inline function source_terms_no_phase_change(u, x, t, equations::CompressibleRainyEulerEquations2D)
    # constants
    g = equations.gravity
    
    # name needed variables
    rho_dry_, rho_moist_, rho_rain_,
    _, rho_v2,
    _, _, _, _,
    rho_dry_h_0, rho_moist_h_0, rho_rain_h_0,
    _, _, _      = u

    # densities
    rho_dry   = rho_dry_    + rho_dry_h_0
    rho_moist = rho_moist_  + rho_moist_h_0
    rho_rain  = rho_rain_   + rho_rain_h_0
    
    rho       = rho_dry     + rho_moist     + rho_rain
    rho_h_0   = rho_dry_h_0 + rho_moist_h_0 + rho_rain_h_0
    rho_      = rho         - rho_h_0

    return SVector(0.0, 0.0, 0.0, 0.0,
                   rho_ * g, 0.0, rho_v2 * g,     #TODO check signs
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
end


@inline function max_abs_speeds(u, equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    v1, v2, v_sound, v_r = cons2speeds(u, equations)

    return abs(v1) + v_sound, abs(v2) + v_sound + abs(v_r)
end


@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleRainyEulerEquations2D)
    # name needed variables
    v1_ll, v2_ll, v_sound_ll, v_r_ll = cons2speeds(u_ll, equations)
    v1_rr, v2_rr, v_sound_rr, v_r_rr = cons2speeds(u_rr, equations)

    # calculate upper bounds for left and right speed
    v_ll_max  = abs(v1_ll * normal_direction[1] +  v2_ll * normal_direction[2])
    v_ll_max += abs(                              v_r_ll * normal_direction[2])
    v_ll_max += v_sound_ll

    v_rr_max  = abs(v1_rr * normal_direction[1] +  v2_rr * normal_direction[2])
    v_rr_max += abs(                              v_r_rr * normal_direction[2])
    v_rr_max += v_sound_rr

    return max(v_ll_max, v_rr_max)
end



###  boundary conditions  ###

@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector, x, t, 
                                              surface_flux_function, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end


@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector, direction, x, t,
                                              surface_flux_function, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end



end  # muladd end