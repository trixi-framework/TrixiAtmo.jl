@muladd begin
#! format: noindent

@doc raw"""
    EnergyTotal
"""
struct EnergyTotal{ThermodynamicStateType, RealType, NGAS} <:
       AbstractThermodynamicEquation
    td_state::ThermodynamicStateType
    inv_gamma_minus_one::SVector{NGAS, RealType}

    function EnergyTotal(td_state)
        @unpack gamma_gas = td_state
        inv_gamma_minus_one = inv.(gamma_gas .- 1)
        return new{typeof(td_state), real(td_state), n_gas(td_state)}(td_state,
                                                                      inv_gamma_minus_one)
    end
end

varname_td(::typeof(cons2cons), ::EnergyTotal) = "rho_e_total"
varname_td(::typeof(cons2prim), ::EnergyTotal) = "p"

# TODO: T_ref!
@inline function temperature(u, equations::CompressibleEulerAtmo{NDIMS}, ::EnergyTotal,
                             td_state::Mixture) where {NDIMS}
    @unpack cv_gas, c_condens, latent_heat = td_state

    # rho_total includes all species
    rho_total = density_total(u, equations)
    rho_gas = vars_gas(u, equations)
    rho_vapor = var_vapor(u, equations)
    rho_condens = vars_condens(u, equations)
    rho_moment = vars_moment(u, equations)

    # Kinetic energy times rho
    rho_Ekin = 0.5f0 * dot(rho_moment, rho_moment) / rho_total

    # Latent energy
    rho_Elat = latent_heat * rho_vapor

    # Absolute temperature
    rho_cv = rho_cv_total(rho_gas, rho_condens, td_state)

    return (var_td(u, equations) - rho_Ekin - rho_Elat) / rho_cv
end

# TODO use other formula to spare division in temperature?
@inline function pressure(u, equations::CompressibleEulerAtmo, ::EnergyTotal,
                          td_state::Mixture)
    @unpack R_gas = td_state

    rho_gas = vars_gas(u, equations)

    T = temperature(u, equations)

    # Gas constant of mixture times rho
    rho_R = rho_R_total(rho_gas, td_state)

    return rho_R * T
end

@inline function speed_of_sound(u, equations::CompressibleEulerAtmo, ::EnergyTotal,
                                td_state::Mixture)
    @unpack R_gas = td_state

    rho_gas = vars_gas(u, equations)
    rho_condens = vars_condens(u, equations)

    cv = cv_total(rho_gas, rho_condens, td_state)
    R = R_total(rho_gas, td_state)
    cp = cv + R
    T = temperature(u, equations)

    return sqrt(cp / cv * R * T)
end

# TODO T_ref
function prim2cons_td(prim, equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS},
                      ::EnergyTotal, td_state::Mixture) where
         {NDIMS, NVARS, NGAS}
    rho_total = prim2density_total(prim, equations)
    prim_gas = vars_gas(prim, equations)
    prim_vapor = var_vapor(prim, equations)
    prim_condens = vars_condens(prim, equations)
    prim_moment = vars_moment(prim, equations)

    cons_gas = SVector{NGAS}(prim_gas[1],
                             (prim_gas[SVector{NGAS - 1}(2:NGAS)] * rho_total)...)

    # calculate temperature
    T = var_td(prim, equations) / rho_R_total(cons_gas, td_state)

    # calculate energy
    Eint = cv_total(cons_gas, prim_condens * rho_total, td_state) * T

    # latent heats
    # TODO fraction only?
    Elatent = td_state.latent_heat * prim_vapor

    # kinetic energy
    Ekinetic = 0.5f0 * dot(prim_moment, prim_moment)

    return (Eint + Elatent + Ekinetic) * rho_total
end

@inline function cons2entropy(cons, equations::CompressibleEulerAtmo{NDIMS},
                              td_equation::EnergyTotal,
                              td_state::Mixture{ParametersType, 1, 0, 0}) where {NDIMS,
                                                                                 ParametersType
                                                                                 }
    rho = density_total(cons, equations)
    p = pressure(cons, equations)
    v = vars_moment(cons, equations) ./ rho

    gamma = td_state.gamma_gas[1]
    s = log(p) - gamma * log(rho)
    rho_p = rho / p

    w1 = (gamma - s) * td_equation.inv_gamma_minus_one[1] -
         0.5f0 * rho_p * dot(v, v)
    w5 = -rho_p

    return SVector(w1, (v * rho_p)..., w5)
end

@inline function flux_td(u, equations::CompressibleEulerAtmo{NDIMS}, ::EnergyTotal,
                         ::Mixture) where {NDIMS}
    p = pressure(u, equations)
    return u[NDIMS + 1] + p
end

@inline function flux_lmars_td(p, v, ::CompressibleEulerAtmo, ::EnergyTotal)
    return p * v
end

@inline function flux_ranocha_td(f_mass_total, rho_total_ll, rho_total_rr, p_ll, p_rr,
                                 v_dot_n_ll, v_dot_n_rr, velocity_square_avg,
                                 ::CompressibleEulerAtmo,
                                 td_equation::EnergyTotal{ThermodynamicStateType,
                                                          RealType, 1}) where {
                                                                               ThermodynamicStateType,
                                                                               RealType}
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_total_ll * p_rr, rho_total_rr * p_ll)

    return f_mass_total *
           (velocity_square_avg + inv_rho_p_mean * td_equation.inv_gamma_minus_one[1]) +
           0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll)
end

@inline function flux_kennedy_gruber_td(f_mass_total, e_avg, p_avg, v_dot_n_avg,
                                        ::CompressibleEulerAtmo, ::EnergyTotal)
    return f_mass_total * e_avg + p_avg * v_dot_n_avg
end

@inline function flux_nonconservative_waruszewski_etal_td(u_ll, u_rr, normal_direction,
                                                          rho_avg, phi_jump,
                                                          equations::CompressibleEulerAtmo,
                                                          ::EnergyTotal)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_avg = 0.5f0 .* (v_ll + v_rr)

    return dot(normal_direction, v_avg) * rho_avg * phi_jump
end

@inline function flux_nonconservative_artiano_etal_td(u_ll, u_rr, normal_direction,
                                                      rho_avg, phi_jump,
                                                      equations::CompressibleEulerAtmo,
                                                      ::EnergyTotal)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_avg = 0.5f0 .* (v_ll + v_rr)

    return dot(normal_direction, v_avg) * rho_avg * phi_jump
end

@inline function flux_nonconservative_souza_etal_td(u_ll, u_rr, normal_direction,
                                                    rho_avg, phi_jump,
                                                    equations::CompressibleEulerAtmo,
                                                    ::EnergyTotal)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_avg = 0.5f0 .* (v_ll + v_rr)

    return dot(normal_direction, v_avg) * rho_avg * phi_jump
end
end # @muladd
