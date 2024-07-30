using Trixi
import Trixi: varnames,
              cons2prim, cons2entropy,
              flux,
              max_abs_speeds, max_abs_speed_naive,
              boundary_condition_slip_wall



#TODO add description of equations and idea + sources



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



###  varnames  ###

varnames(::typeof(cons2cons), ::CompressibleRainyEulerEquations2D) = ("_rho_dry", "_rho_moist", "_rho_rain",
                                                                      "rho_v1", "rho_v2",
                                                                      "_energy",
                                                                      "rho_vapour", "rho_cloud",
                                                                      "temperature",
                                                                      "rho_dry_hydrostatic", "rho_moist_hydrostatic", "rho_rain_hydrostatic",
                                                                      "energy_hydrostatic",
                                                                      "rho_vapour_hydrostatic",
                                                                      "temperature_hydrostatic")


varnames(::typeof(cons2prim), ::CompressibleRainyEulerEquations2D) = ("rho_dry", "rho_moist", "rho_rain",
                                                                      "v1", "v2",
                                                                      "energy",
                                                                      "rho_vapour", "rho_cloud",
                                                                      "temperature",
                                                                      "rho_dry_hydrostatic", "rho_moist_hydrostatic", "rho_rain_hydrostatic",
                                                                      "energy_hydrostatic",
                                                                      "rho_vapour_hydrostatic",
                                                                      "temperature_hydrostatic")




###  conversion  ###

@inline function cons2prim(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return SVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end


@inline function cons2entropy(u, equations::CompressibleRainyEulerEquations2D)
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


@inline function flux(u, normal_direction::AbstractVector, equations::CompressibleRainyEulerEquations2D)
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


@inline function max_abs_speeds(u, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return 0.0, 0.0
end


@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleRainyEulerEquations2D)
    #TODO
    return 0.0
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