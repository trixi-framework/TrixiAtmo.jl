# Disable formatting, keep this file nicely looking!
#! format: off

abstract type AbstractParameters end

@doc raw"""
    Parameters for TrixiAtmo simulations

# TODO
       
- some things depend on temperature, make functions instead?
- aliases to have shorter names?

"""
@kwdef struct Parameters{RealType <: Real} <: AbstractParameters

    # Specific heat capacities
    c_dry_air_const_pressure::RealType         = RealType(1004.5)    # Specific heat of dry air at constant pressure
    c_dry_air_const_volume::RealType           = RealType(717.5)     # Specific heat of dry air at constant volume
    c_vapour_const_pressure:: RealType         = RealType(1859.5)    # Specific heat of vapour at constant pressure
    c_vapour_const_volume:: RealType           = RealType(1397.5)    # Specific heat of vapour at constant volume
    c_liquid_water:: RealType                  = RealType(4181.0)    # Specific heat of water

    # Reference values
    ref_pressure::RealType                     = RealType(1e5)
    ref_temperature::RealType                  = RealType(273.16)
    ref_latent_heat_vaporization::RealType     = RealType(2.5008e6)  # for ref_temperature!
    ref_saturation_pressure::RealType          = RealType(611.657)   # for ref_temperature!

    # Earth
    earth_radius::RealType                     = RealType(6.37122e6) # m
    earth_gravitational_acceleration::RealType = RealType(9.80616)   # m/s²
    earth_rotation_rate::RealType              = RealType(7.292e-5)  # rad/s

    # Other:
    seconds_per_day::RealType                  = RealType(8.64e4)

    # Numerics
    tol_eps::RealType                          = eps(RealType)       # used to prevent zeros in some calculations, can be overridden to disable
end
