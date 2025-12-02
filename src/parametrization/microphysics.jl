@muladd begin
#! format: noindent

abstract type AbstractMicroPhysics end

@doc raw"""
    MicrophysicsRelaxation( ??? )

"""
struct MicrophysicsRelaxation{RealType <: Real} <: AbstractMicroPhysics
    velocity_rain_constant::RealType

    function MicrophysicsRelaxation{RealType}() where {RealType}
        # COSMO model
        cosmo_constant = (pi * 8f6)^(-0.125) * gamma(4.5f0) * 130 / 6
        return new{RealType}(cosmo_constant)
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
end # @muladd

