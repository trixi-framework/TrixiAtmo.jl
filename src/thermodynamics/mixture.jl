@muladd begin
#! format: noindent

############################################################################################
# AbstractThermodynamicState
############################################################################################

abstract type AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP} end

varnames_gas(::Union{typeof(cons2cons), typeof(cons2prim)},
::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}) where
{RealType, NGAS, NCONDENS, NPRECIP} = ntuple(i -> "gas_$i", Val(NGAS))

varnames_liquid(::Union{typeof(cons2cons), typeof(cons2prim)},
::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}) where
{RealType, NGAS, NCONDENS, NPRECIP} = ntuple(i -> "liquid_$i", Val(NCONDENS))

varnames_precip(::Union{typeof(cons2cons), typeof(cons2prim)},
::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}) where
{RealType, NGAS, NCONDENS, NPRECIP} = ntuple(i -> "precip_$i", Val(NPRECIP))

@inline function Base.real(::AbstractThermodynamicState{RealType}) where {RealType}
    return RealType
end

############################################################################################
# Mixture
############################################################################################

@doc raw"""
    Mixture
"""
struct Mixture{ParametersType, NGAS, NCONDENS, NPRECIP, RealType} <:
       AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP}

    #gas_names::GasNames
    #condens_names::CondensNames
    #precip_names::PrecipNames

    parameters::ParametersType
    cv_gas::SVector{NGAS, RealType}
    cp_gas::SVector{NGAS, RealType}
    R_gas::SVector{NGAS, RealType}
    gamma_gas::SVector{NGAS, RealType}
    c_condens::SVector{NCONDENS, RealType}

    latent_heat::RealType # TODO size?
    p_ref::RealType

    varnames_gas::SVector{NGAS, String}
    varnames_liquid::SVector{NCONDENS, String}
    varnames_precip::SVector{NPRECIP, String}

    function Mixture(;
                     parameters::Parameters{RealType},
                     gases::Tuple = (DryAir(parameters),),
                     condensates::Tuple = (),
                     precipitates::Tuple = ()) where {RealType}
        n_gas = length(gases)
        n_condens = length(condensates)
        n_precip = length(precipitates)

        cv_gas = SVector{n_gas}(gas.cv for gas in gases)
        cp_gas = SVector{n_gas}(gas.cp for gas in gases)
        R_gas = SVector{n_gas}(gas.R for gas in gases)
        gamma_gas = SVector{n_gas}(gas.gamma for gas in gases)

        c_condens = SVector{n_condens}(condens.c for condens in condensates)

        # TODO ?
        latent_heat = parameters.ref_latent_heat_vaporization

        p_ref = parameters.ref_pressure

        varnames_gas = SVector{n_gas}(gas.varname for gas in gases)
        varnames_condens = SVector{n_condens}(condens.varname
                                              for condens in condensates)
        varnames_precip = SVector{n_precip}(precip.varname for precip in precipitates)

        return new{typeof(parameters), n_gas, n_condens, n_precip, RealType}(parameters,
                                                                             cv_gas,
                                                                             cp_gas,
                                                                             R_gas,
                                                                             gamma_gas,
                                                                             c_condens,
                                                                             latent_heat,
                                                                             p_ref,
                                                                             varnames_gas,
                                                                             varnames_condens,
                                                                             varnames_precip)
    end
end

varnames_gas(::Union{typeof(cons2cons), typeof(cons2prim)},
td_state::Mixture) = td_state.varnames_gas

varnames_liquid(::Union{typeof(cons2cons), typeof(cons2prim)},
td_state::Mixture) = td_state.varnames_liquid

varnames_precip(::Union{typeof(cons2cons), typeof(cons2prim)},
td_state::Mixture) = td_state.varnames_precip

n_gas(::Mixture{ParametersType, NGAS}) where {ParametersType, NGAS} = NGAS

# used in equation of state: p = rho R T    
# pressure is sum of partial pressures, usually of dry air and water vapor
# neglects condensed phases (nearly no volume, but mass)
@inline function rho_R_total(rho_gas, td_state::Mixture)
    @unpack R_gas = td_state
    return dot(rho_gas, td_state.R_gas)
end
@inline function R_total(rho_gas, td_state::Mixture)
    rho_total = sum(rho_gas)
    return rho_R_total(rho_gas, td_state) / rho_total
end

# used in internal energy
# TODO: this includes condensates
@inline function rho_cv_total(rho_gas, rho_condens, td_state::Mixture)
    @unpack cv_gas, c_condens = td_state
    return dot(rho_gas, cv_gas) + dot(rho_condens, c_condens)
end
@inline function cv_total(rho_gas, rho_condens, td_state::Mixture)
    rho_total = sum(rho_gas) + sum(rho_condens)
    return rho_cv_total(rho_gas, rho_condens, td_state) / rho_total
end

@inline function gamma_total(rho_gas, rho_condens, td_state::Mixture)
    error("gamma_total missing multi species implementation")
end

#=
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
=#

############################################################################################
# IdealGas (special case)
############################################################################################

const IdealGas{RealType} = Mixture{Parameters, 1, 0, 0, RealType}
IdealGas(; parameters) = Mixture(; parameters)

@inline function gamma_total(rho_gas, rho_condens,
                             td_state::Mixture{ParametersType, 1, 0, 0}) where {ParametersType}
    return td_state.gamma_gas[1]
end
end # @muladd
