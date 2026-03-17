@muladd begin
#! format: noindent

abstract type AbstractMicroPhysics end

@doc raw"""
    MicrophysicsRelaxation( ??? )

"""
struct MicrophysicsRelaxation{RealType <: Real} <: AbstractMicroPhysics
    saturation_factor::RealType
    velocity_rain_constant::RealType

    function MicrophysicsRelaxation{RealType}(;
        saturation_factor = RealType(1)) where {RealType}
        # COSMO model
        cosmo_constant = (pi * 8f6)^(-0.125) * gamma(4.5f0) * 130 / 6
        return new{RealType}(saturation_factor, cosmo_constant)
    end
end

@inline function velocity_rain(u, equations, 
                               microphysics::MicrophysicsRelaxation)
    @unpack velocity_rain_constant = microphysics
    rho_airborn = vars_airborn(u, equations)
    rho_precip = vars_precip(u, equations)

    # TODO: implicit assumption: rain comes first
    rho_rain = rho_precip[1]

    # else problems when rho_water == 0 as well
    if rho_rain == 0
        return rho_rain
    end

    # TODO: all airborn minus dry air, is this even correct?
    rho_water = sum(rho_airborn) - rho_airborn[1]

    return velocity_rain_constant *
           (rho_rain / rho_water)^(0.125f0)
end

@inline function saturation_vapour_pressure(temperature,
                                            parameters,
                                            ::MicrophysicsRelaxation)
    c_l = parameters.c_liquid_water
    c_pv = parameters.c_vapour_const_pressure
    c_vv = parameters.c_vapour_const_volume
    R_v = c_pv - c_vv

    ref_s_p = parameters.ref_saturation_pressure
    ref_temp = parameters.ref_temperature
    ref_L = parameters.ref_latent_heat_vaporization

    # Clausius Clapeyron formula
    return ref_s_p * (temperature / ref_temp)^((c_pv - c_l) / R_v) *
           exp(((ref_L - (c_pv - c_l) * ref_temp) / R_v) * (1 / ref_temp - 1 / temperature))
end

# This source term models the phase chance between could water and vapor.
@inline function phase_change_vapour_liquid(u, equations,
                                            microphysics::MicrophysicsRelaxation)
    @unpack saturation_factor = microphysics

    c_pv = equations.parameters.c_vapour_const_pressure
    c_vv = equations.parameters.c_vapour_const_volume
    R_v = c_pv - c_vv
    
    rho_d, rho_v = vars_gas(u, equations)
    rho = density_total(u, equations)
    T = temperature(u, equations)

    # saturation vapor pressure
    p_vs = saturation_vapour_pressure(T, equations.parameters, microphysics)

    # saturation density of vapor 
    rho_star_qv = p_vs / (R_v * T)

    # Fisher-Burgmeister-Function
    a = rho_star_qv - rho_v
    b = rho - rho_v - rho_d

    return (a + b - sqrt(a^2 + b^2)) * saturation_factor
end
end # @muladd
