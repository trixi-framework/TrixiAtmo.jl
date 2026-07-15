@muladd begin
#! format: noindent

@doc raw"""
    PotentialTemperature
"""
struct PotentialTemperature{ThermodynamicStateType, RealType, NGAS} <:
       AbstractThermodynamicEquation
    td_state::ThermodynamicStateType
    K::SVector{NGAS, RealType} # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    #stolarsky_factor::SVector{NGAS, RealType} # = (gamma - 1) / gamma; used in the stolarsky mean

    function PotentialTemperature(td_state)
        @unpack R_gas, gamma_gas, p_ref = td_state
        K = p_ref * (R_gas / p_ref) .^ gamma_gas
        #stolarsky_factor = (gamma_gas .- 1) / gamma_gas
        return new{typeof(td_state), real(td_state), n_gas(td_state)}(td_state, K)
        #stolarsky_factor)
    end
end

varname_td(::typeof(cons2cons), ::PotentialTemperature) = "rho_theta"
varname_td(::typeof(cons2prim), ::PotentialTemperature) = "p"

@inline function pressure(u, equations::CompressibleEulerAtmo{NDIMS},
                          td_equation::PotentialTemperature,
                          td_state::Mixture) where {NDIMS}
    gas = vars_gas(u, equations)
    condens = vars_condens(u, equations)
    K = K_total(gas, condens, td_equation, td_state)
    gamma = gamma_total(gas, condens, td_state)
    # way faster than using ^td_state.gamma_gas[1]
    p = K * exp(gamma * log(var_td(u, equations)))
    return p
end

@inline function speed_of_sound(u, equations::CompressibleEulerAtmo{NDIMS},
                                td_equation::PotentialTemperature,
                                td_state::Mixture) where {NDIMS}
    p = pressure(u, equations, td_equation, td_state)
    gamma = gamma_total(vars_gas(u, equations), vars_condens(u, equations), td_state)
    return sqrt(gamma * p / u[NDIMS + 2])
end

@inline function K_total(rho_gas, rho_condens, td_equation::PotentialTemperature,
                         td_state::Mixture{ParametersType, 1, 0, 0}) where {ParametersType}
    return td_equation.K[1]
end

function prim2cons_td(prim, equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS},
                      td_equation::PotentialTemperature,
                      td_state::Mixture) where {NDIMS, NVARS, NGAS}
    @unpack p_ref = td_state

    rho_total = prim2density_total(prim, equations)
    prim_gas = vars_gas(prim, equations)
    condens = vars_condens(prim, equations)
    cons_gas = SVector{NGAS}(prim_gas[1],
                             (prim_gas[SVector{NGAS - 1}(2:NGAS)] * rho_total)...)
    K = K_total(cons_gas, condens, td_equation, td_state)
    gamma = gamma_total(cons_gas, condens * rho_total, td_state)

    rho_theta = (var_td(prim, equations) / K)^(1 / gamma)
    return rho_theta
end

@inline function cons2entropy(cons, ::CompressibleEulerAtmo{NDIMS, NVARS},
                              td_equation::PotentialTemperature,
                              td_state::Mixture{ParametersType, 1, 0, 0}) where {NDIMS, NVARS, ParametersType}
    rho = density_total(cons, equations)
    rho_theta = var_td(cons, equations)

    K = td_equation.K[1]
    gamma = td_state.gamma_gas[1]

    w1 = log(K * (rho_theta / rho)^gamma) - gamma
    w5 = rho / rho_theta * gamma

    n_passive = NVARS - NDIMS - 2

    return SVector(ntuple(i -> 0, Val(NDIMS))..., w5, w1, ntuple(i -> 0, Val(n_passive))...)
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

# Only implemented for dry air. Mixture could results in gamma_ll and gamma_rr.
@inline function flux_tec_td(rho_theta_ll, rho_theta_rr, v_dot_n_avg,
                             equations::CompressibleEulerAtmo,
                             ::PotentialTemperature{ThermodynamicStateType, RealType,
                                                    1}) where {ThermodynamicStateType,
                                                               RealType}
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr,
                               equations.td_state.gamma_gas[1])
    return gammamean * v_dot_n_avg
end

# Only implemented for dry air. Mixture could results in gamma_ll and gamma_rr.
@inline function flux_etec_td(rho_theta_ll, rho_theta_rr, v_dot_n_avg,
                              equations::CompressibleEulerAtmo,
                              ::PotentialTemperature{ThermodynamicStateType, RealType,
                                                     1}) where {ThermodynamicStateType,
                                                                RealType}
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr,
                               equations.td_state.gamma_gas[1])
    return gammamean * v_dot_n_avg
end

@inline function flux_nonconservative_waruszewski_etal_td(u_ll, u_rr, normal_direction,
                                                          rho_avg, phi_jump,
                                                          ::CompressibleEulerAtmo,
                                                          ::PotentialTemperature)
    return 0
end

@inline function flux_nonconservative_artiano_etal_td(u_ll, u_rr, normal_direction,
                                                      rho_avg, phi_jump,
                                                      ::CompressibleEulerAtmo,
                                                      ::PotentialTemperature)
    return 0
end

@inline function flux_nonconservative_souza_etal_td(u_ll, u_rr, normal_direction,
                                                    rho_avg, phi_jump,
                                                    ::CompressibleEulerAtmo,
                                                    ::PotentialTemperature)
    return 0
end
end # @muladd
