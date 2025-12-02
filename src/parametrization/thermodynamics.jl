@muladd begin
#! format: noindent

############################################################################################
# AbstractThermodynamicState
############################################################################################

abstract type AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP} end

varnames_gas(::Union{typeof(cons2cons), typeof(cons2prim)},
             ::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}) where
             {RealType, NGAS, NCONDENS, NPRECIP} =
    ntuple(i -> "_gas_$i", Val(NGAS))

varnames_liquid(::Union{typeof(cons2cons), typeof(cons2prim)},
                ::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}) where
                {RealType, NGAS, NCONDENS, NPRECIP} =
    ntuple(i -> "_liquid_$i", Val(NCONDENS))

varnames_precip(::Union{typeof(cons2cons), typeof(cons2prim)},
                ::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}) where
                {RealType, NGAS, NCONDENS, NPRECIP} =
    ntuple(i -> "_precip_$i", Val(NPRECIP))

@inline function Base.real(::AbstractThermodynamicState{RealType}) where {RealType}
    return RealType
end

############################################################################################
# IdealGas
############################################################################################

struct IdealGas{RealType} <: AbstractThermodynamicState{RealType, 1, 0, 0}
    cv::RealType
    cp::RealType
    gas_constant::RealType
    gamma::RealType
    p_ref::RealType

    function IdealGas(; parameters::Parameters{RealType}) where {RealType}
        cv = parameters.c_dry_air_const_volume
        cp = parameters.c_dry_air_const_pressure
        R = cp - cv
        gamma = cp / cv
        p_ref = parameters.ref_pressure
        return new{RealType}(cv, cp, R, gamma, p_ref)
    end
end

@inline function gamma_total(rho_gas, rho_condens, td_state::IdealGas)
    return td_state.gamma
end

varnames_gas(::typeof(cons2cons), ::IdealGas) = ("", )



############################################################################################
# IdealGasesAndLiquids
############################################################################################

@doc raw"""

This designed to be in line with CompressibleEulerAtmo

Select how many of the following predefined will be used:

Gaseous species (n_gas)
1. dry air
2. vapor 

Condensed species (n_condens)
1. cloud water

Precipitating species (n_precip)
1. rain

"""
struct IdealGasesAndLiquids{ParametersType, NGAS, NCONDENS, NPRECIP, RealType} <:
    AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}

    parameters::ParametersType
    cv_gas::SVector{NGAS, RealType}
    cp_gas::SVector{NGAS, RealType}
    gas_constant::SVector{NGAS, RealType}
    c_condens::SVector{NCONDENS, RealType}
    latent_heat::SVector{NCONDENS, RealType}

    varnames_gas::SVector{NGAS, String}
    varnames_liquid::SVector{NCONDENS, String}
    varnames_precip::SVector{NPRECIP, String}

    function IdealGasesAndLiquids(;
        parameters::Parameters{RealType},
        n_gas = 1, n_condens = 0, n_precip = 0) where {RealType}

        @assert n_gas <=2 "Only up to 2 gaseous species supported"
        @assert n_condens <=1 "Only up to 1 condensed species supported"
        @assert n_precip <=1  "Only up to 1 precipitating species supported"

        cv_gas_params = [:c_dry_air_const_volume, :c_vapour_const_volume]
        cp_gas_params = [:c_dry_air_const_pressure, :c_vapour_const_pressure]
        c_condens_param = [:c_liquid_water]
        latent_heat = [:latent_heat_vaporization]

        cv_gas = SVector{n_gas}(getproperty(parameters, cv_gas_params[i]) for i in 1:n_gas)
        cp_gas = SVector{n_gas}(getproperty(parameters, cp_gas_params[i]) for i in 1:n_gas)
        c_condens = SVector{n_condens}(getproperty(parameters, c_condens_param[i])
                                        for i in 1:n_condens)
        latent_heat = SVector{n_condens}(getproperty(parameters, latent_heat[i])
                                        for i in 1:n_condens)
        gas_constant = cp_gas .- cv_gas

        varnames_gas = SVector{n_gas}(["_dry", "_vapor"][1:n_gas])
        varnames_condens = SVector{n_condens}(["_cloud"][1:n_condens])
        varnames_precip = SVector{n_precip}(["_rain"][1:n_precip])

        return new{typeof(parameters), n_gas, n_condens, n_precip, RealType}(
            parameters, cv_gas, cp_gas, gas_constant, c_condens, latent_heat,
            varnames_gas, varnames_condens, varnames_precip)
    end
end

varnames_gas(::typeof(cons2cons),
             td_state::IdealGasesAndLiquids) = td_state.varnames_gas

varnames_liquid(::typeof(cons2cons),
                td_state::IdealGasesAndLiquids) = td_state.varnames_liquid

varnames_precip(::typeof(cons2cons),
                td_state::IdealGasesAndLiquids) = td_state.varnames_precip

# used in equation of state: p = rho R T    
# pressure is sum of partial pressures, usually of dry air and water vapor
# neglects condensed phases (nearly no volume, but mass)
@inline function R_total(rho_gas, td_state::IdealGasesAndLiquids)
    rho_total = sum(rho_gas)
    return dot(rho_gas, td_state.gas_constant) / rho_total
end

# used in internal energy
@inline function cv_total(rho_gas, rho_condens, td_state::IdealGasesAndLiquids)
    @unpack cv_gas, c_condens = td_state
    rho_total = sum(rho_gas) + sum(rho_condens)
    return (dot(rho_gas, cv_gas) + dot(rho_condens, c_condens)) / rho_total
end

@inline function gamma_total(rho_gas, rho_condens,
    td_state::IdealGasesAndLiquids{RealType, NGAS}) where {RealType, NGAS}

    return R_total(rho_gas, td_state) /
           cv_total(rho_gas, rho_condens, td_state) + 1
end

@inline function entropies_gas(densities, T,
    td_state::IdealGasesAndLiquids{RealType, NGAS}) where {RealType, NGAS}

    @unpack cp_gas, gas_constant, parameters = td_state

    # TODO -> ctor
    s_ref = SVector{NGAS}(0 for i in 1:NGAS)

    T_ref = parameters.ref_temperature
    p_ref = parameters.ref_pressure
    eps = parameters.tol_eps
    
    # partial pressures
    p = densities .* gas_constant .* T .+ eps

    return s_ref .+ cp_gas .* log(T / T_ref) .- gas_constant .* log.(p / p_ref)
end

@inline function entropies_liquid(densities, T,
    td_state::IdealGasesAndLiquids{RealType, NGAS, NCONDENS}) where {RealType, NGAS, NCONDENS}

    @unpack parameters = td_state

    # TODO -> ctor, other latent
    latent = SVector{NCONDENS}(parameters.latent_heat_vaporization for i in 1:NCONDENS)

    return latent ./ T
end
end # @muladd
