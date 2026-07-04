@muladd begin
#! format: noindent

@doc raw"""
    PotentialTemperature
"""
struct PotentialTemperature{ThermodynamicStateType, RealType, NGAS} <:
       AbstractThermodynamicEquation
    td_state::ThermodynamicStateType
    K::SVector{NGAS, RealType} # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::SVector{NGAS, RealType} # = (gamma - 1) / gamma; used in the stolarsky mean

    function PotentialTemperature(td_state)
        @unpack R_gas, gamma_gas, p_ref = td_state
        K = p_ref * (R_gas / p_ref) .^ gamma_gas
        stolarsky_factor = (gamma_gas .- 1) / gamma_gas
        return new{typeof(td_state), real(td_state), n_gas(td_state)}(td_state, K,
                                                                      stolarsky_factor)
    end
end

varname_td(::typeof(cons2cons), ::PotentialTemperature) = "rho_theta"
varname_td(::typeof(cons2prim), ::PotentialTemperature) = "p"

@inline function pressure(u, equations::CompressibleEulerAtmo{NDIMS},
                          td_equation::PotentialTemperature,
                          td_state::Mixture) where {NDIMS}
    p = td_equation.K[1] * var_td(u, equations)^td_state.gamma_gas[1]
    return p
end

@inline function speed_of_sound(u, equations::CompressibleEulerAtmo{NDIMS},
                                td_equation::PotentialTemperature,
                                td_state::Mixture) where {NDIMS}
    p = pressure(u, equations, td_equation, td_state)
    # TODO: single only, total density?
    return sqrt(td_state.gamma_gas[1] * p / u[NDIMS + 2])
end

@inline function flux_td(u, equations::CompressibleEulerAtmo{NDIMS},
                         ::PotentialTemperature,
                         ::Mixture) where {NDIMS}
    return var_td(u, equations)
end

@inline function flux_lmars_td(p, v, ::CompressibleEulerAtmo{NDIMS, NVARS},
                               ::PotentialTemperature) where {NDIMS, NVARS}
    return 0
end

@inline function flux_ec_td(f_mass_total, rho_theta_ll, rho_theta_rr, rho_total_ll,
                            rho_total_rr,
                            equations::CompressibleEulerAtmo, ::PotentialTemperature)
    return f_mass_total *
           inv_ln_mean(rho_total_ll / rho_theta_ll, rho_total_rr / rho_theta_rr)
end

@inline function flux_tec_td(rho_theta_ll, rho_theta_rr, v_dot_n_avg,
                             equations::CompressibleEulerAtmo, ::PotentialTemperature)
    # TODO only 1
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr,
                               equations.td_state.gamma_gas[1])
    return gammamean * v_dot_n_avg
end

@inline function flux_etec_td(rho_theta_ll, rho_theta_rr, v_dot_n_avg,
                              equations::CompressibleEulerAtmo, ::PotentialTemperature)
    # TODO only 1
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr,
                               equations.td_state.gamma_gas[1])
    return gammamean * v_dot_n_avg
end

function prim2cons_td(prim, equations::CompressibleEulerAtmo{NDIMS},
                      ::PotentialTemperature,
                      td_state::Mixture) where {NDIMS}
    @unpack R_gas, gamma_gas, p_ref = td_state
    # TODO: single only

    return (var_td(prim, equations) / p_ref)^(1 / gamma_gas[1]) * p_ref / R_gas[1]
end
end # @muladd
