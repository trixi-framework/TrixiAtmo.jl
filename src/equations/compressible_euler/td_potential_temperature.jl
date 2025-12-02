@muladd begin
#! format: noindent

@doc raw"""
    PotentialTemperature
"""
struct PotentialTemperature{ThermodynamicStateType, RealType} <: AbstractThermodynamicEquation
    td_state::ThermodynamicStateType
    K::RealType # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealType # = (gamma - 1) / gamma; used in the stolarsky mean

    function PotentialTemperature(td_state)
        @unpack gas_constant, gamma, p_ref = td_state
        K = p_ref * (gas_constant / p_ref)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        return new{typeof(td_state), real(td_state)}(td_state, K, stolarsky_factor)
    end
end

varname_td(::typeof(cons2cons), ::PotentialTemperature) = "rho_theta"
varname_td(::typeof(cons2prim), ::PotentialTemperature) = "p"

@inline function pressure(u, equations::CompressibleEulerAtmo{NDIMS},
                          td_equation::PotentialTemperature,
                          td_state::IdealGas) where {NDIMS}
    p = td_equation.K * exp(td_state.gamma * log(u[NDIMS+1]))
    return p
end

@inline function speed_of_sound(u, equations::CompressibleEulerAtmo{NDIMS},
                                td_equation::PotentialTemperature,
                                td_state::IdealGas) where {NDIMS}
    p = pressure(u, equations, td_equation, td_state)
    return sqrt(td_state.gamma * p / u[NDIMS+2])
end

@inline function flux_td(u, equations::CompressibleEulerAtmo{NDIMS}, ::PotentialTemperature,
                         ::IdealGas) where {NDIMS}
    return u[NDIMS+1]
end

function prim2cons_td(prim, equations::CompressibleEulerAtmo{NDIMS},
                      ::PotentialTemperature,
                      td_state::IdealGas) where {NDIMS}
    @unpack gas_constant, gamma, p_ref = td_state
    return (prim[NDIMS+1] / p_ref)^(1 / gamma) * p_ref / gas_constant
end

end # @muladd
