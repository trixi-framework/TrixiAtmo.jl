@muladd begin
#! format: noindent


############################################################################################
# AbstractSpecies
############################################################################################

abstract type AbstractSpecies{RealType} end
abstract type AbstractGas{RealType} <: AbstractSpecies{RealType} end
abstract type AbstractLiquid{RealType} <: AbstractSpecies{RealType} end
abstract type AbstractPrecipitation{RealType} <: AbstractSpecies{RealType} end


############################################################################################
# Gases
############################################################################################

@doc raw"""
    DryAir
"""
struct DryAir{RealType} <: AbstractGas{RealType}
    cv::RealType
    cp::RealType
    R::RealType
    gamma::RealType
    varname::String

    function DryAir(parameters::Parameters{RealType}) where {RealType}
        cv = parameters.c_dry_air_const_volume
        cp = parameters.c_dry_air_const_pressure
        R = cp - cv
        gamma = cp / cv
        return new{RealType}(cv, cp, R, gamma, "dry")
    end
end

@doc raw"""
    WaterVapor
"""
struct WaterVapor{RealType} <: AbstractGas{RealType}
    cv::RealType
    cp::RealType
    R::RealType
    gamma::RealType
    varname::String

    function WaterVapor(parameters::Parameters{RealType}) where {RealType}
        cv = parameters.c_vapor_const_volume
        cp = parameters.c_vapor_const_pressure
        R = cp - cv
        gamma = cp / cv
        return new{RealType}(cv, cp, R, gamma, "vapor")
    end
end


############################################################################################
# Condensates
############################################################################################

@doc raw"""
    CloudWater
"""
struct CloudWater{RealType} <: AbstractLiquid{RealType}
    c::RealType
    varname::String

    function CloudWater(parameters::Parameters{RealType}) where {RealType}
        c = parameters.c_liquid_water
        return new{RealType}(c, "cloud")
    end
end


############################################################################################
# Precipitates
############################################################################################

@doc raw"""
    RainWater
"""
struct RainWater{RealType} <: AbstractPrecipitation{RealType}
    c::RealType
    varname::String

    function RainWater(parameters::Parameters{RealType}) where {RealType}
        c = parameters.c_liquid_water
        return new{RealType}(c, "rain")
    end
end

end # @muladd