using Trixi



###  equation, parameters and constants  ###

struct CompressibleRainyEulerEquations2D{RealT <: Real}
    # Specific heat capacities:
    c_liquid_water             ::Float64
    c_dry_air_const_pressure   ::Float64
    c_dry_air_const_volume     ::Float64
    c_vapour_const_pressure    ::Float64
    c_vapour_const_volume      ::Float64

    # Gas constants:
    R_dry_air                  ::Float64
    R_vapour                   ::Float64
    eps                        ::Float64

    # Reference values:
    ref_saturation_pressure    ::Float64
    ref_temperature            ::Float64
    ref_latent_heat_vap_temp   ::Float64
    ref_pressure               ::Float64

    # Other:
    gravity                    ::Float64
    rain_water_distr           ::Float64
    v_mean_rain                ::Float64
end


function CompressibleRainyEulerEquations2D()
    # Specific heat capacities:
    c_liquid_water = 4186.0
    c_dry_air_const_pressure = 1004.0
    c_dry_air_const_volume = 717.0
    c_vapour_const_pressure = 1885.0
    c_vapour_const_volume = 1424.0

    # Gas constants:
    R_dry_air = c_dry_air_const_pressure - c_dry_air_const_volume
    R_vapour = c_vapour_const_pressure - c_vapour_const_volume
    eps = R_dry_air / R_vapour

    # Reference values:
    saturation_vapour_pressure = 610.7
    ref_temperature = 273.15
    latent_heat_vap_ref_temp = 2.5e6
    ref_pressure  = 1e5

    # Other:
    gravity = 9.81
    rain_water_distr = 8e6
    v_mean_rain = 130.0

    return CompressibleRainyEulerEquations2D{RealT}(c_liquid_water, c_dry_air_const_pressure, c_dry_air_const_volume, 
                          c_vapour_const_pressure, c_vapour_const_volume, R_dry_air, 
                          R_vapour, eps, saturation_vapour_pressure, ref_temperature,
                          latent_heat_vap_ref_temp, ref_pressure, gravity, 
                          rain_water_distr, v_mean_rain)
end



###  conversion  ###

@inline function cons2prim(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end



###  physics variables  ###

@inline function speed_of_sound(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return v_sound
end

@inline function terminal_velocity_rain(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return v_terminal_rain
end

@inline function saturation_vapour_pressure(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return p_vapour_saturation
end



###  pde discretization  ###

@inline function flux(u, orientation::Integer, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

@inline function source_terms_rainy(u, x, t, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

@inline function source_terms_no_rain(u, x, t, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end