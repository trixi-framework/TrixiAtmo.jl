@muladd begin
#! format: noindent

@doc raw"""
    EnergyInternal
"""
struct EnergyInternal{ThermodynamicStateType} <: AbstractThermodynamicEquation
    td_state::ThermodynamicStateType
end

varname_td(::typeof(cons2cons), ::EnergyInternal) = "rho_e_internal"
varname_td(::typeof(cons2prim), ::EnergyInternal) = "p"

# TODO: T_ref!
@inline function temperature(u, equations::CompressibleEulerAtmo{NDIMS},
                             ::EnergyInternal,
                             td_state::Mixture) where {NDIMS}
    @unpack cv_gas, c_condens, latent_heat = td_state

    rho_total = density_total(u, equations)
    rho_gas = vars_gas(u, equations)
    rho_vapor = var_vapor(u, equations)
    rho_condens = vars_condens(u, equations)

    # Latent energy
    # TODO: ?
    rho_Elat = latent_heat * rho_vapor

    # Absolute temperature
    rho_cv = rho_total * cv_total(rho_gas, rho_condens, td_state)

    return (var_td(u, equations) - rho_Elat) / rho_cv
end

@inline function pressure(u, equations::CompressibleEulerAtmo, ::EnergyInternal,
                          td_state::Mixture)
    rho_gas = vars_gas(u, equations)
    rho_condens = vars_condens(u, equations)

    # TODO
    #rho_Elat = latent_heat * rho_vapor

    gamma = gamma_total(rho_gas, rho_condens, td_state)

    return (gamma - 1) * var_td(u, equations)
end

@inline function speed_of_sound(u, equations::CompressibleEulerAtmo, ::EnergyInternal,
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
                      ::EnergyInternal, td_state::Mixture) where
         {NDIMS, NVARS, NGAS}
    rho_total = prim2density_total(prim, equations)
    prim_gas = vars_gas(prim, equations)
    prim_vapor = var_vapor(prim, equations)
    prim_condens = vars_condens(prim, equations)

    cons_gas = SVector{NGAS}(prim_gas[1],
                             (prim_gas[SVector{NGAS - 1}(2:NGAS)] * rho_total)...)

    gamma = gamma_total(cons_gas, prim_condens * rho_total, td_state)

    # calculate energy
    Eint = var_td(prim, equations) / (gamma - 1)

    # latent heats
    # TODO fraction only?
    Elatent = td_state.latent_heat * prim_vapor

    return Eint + Elatent
end

@inline function cons2entropy(cons, equations::CompressibleEulerAtmo{NDIMS, NVARS},
                              ::EnergyInternal,
                              td_state::Mixture{ParametersType, 1, 0, 0}) where {NDIMS, NVARS, ParametersType}
    rho = density_total(cons, equations)
    rho_e = var_td(cons, equations)

    gamma = td_state.gamma_gas[1]

    w1 = log(rho_e * (gamma - 1) / rho^gamma) - gamma
    w5 = rho / (rho_e * (gamma - 1))

    n_passive = NVARS - NDIMS - 2

    return SVector(zeros(SVector{3, Int64})..., w5, w1, ntuple(i -> 0, Val(n_passive))...)
end

@inline function flux_td(u, equations::CompressibleEulerAtmo{NDIMS}, ::EnergyInternal,
                         ::Mixture) where {NDIMS}
    p = pressure(u, equations)
    return u[NDIMS + 1] + p
end

@inline function flux_lmars_td(p, v, ::CompressibleEulerAtmo, ::EnergyInternal)
    return p * v
end

@inline function flux_kennedy_gruber_td(f_mass_total, e_avg, p_avg, v_dot_n_avg,
                                        ::CompressibleEulerAtmo, ::EnergyInternal)
    return f_mass_total * e_avg + p_avg * v_dot_n_avg
end

@inline function flux_artiano_td(rho_ll, rho_rr, p_ll, p_rr,
                                 f_mass_total, v_interface,
                                 equations::CompressibleEulerAtmo,
                                 td_equation::EnergyInternal)
    # TODO only one!
    T_ll = p_ll / (td_equation.td_state.R_gas[1] * rho_ll)
    T_rr = p_rr / (td_equation.td_state.R_gas[1] * rho_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    T_avg = inv_ln_mean(1 / T_ll, 1 / T_rr) * td_equation.td_state.cv_gas[1]

    return f_mass_total * T_avg + p_avg * v_interface
end

@inline function flux_nonconservative_waruszewski_etal_td(u_ll, u_rr, normal_direction,
                                                          rho_avg, phi_jump,
                                                          equations::CompressibleEulerAtmo,
                                                          ::EnergyInternal)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_avg = 0.5f0 .* (v_ll + v_rr)

    return -dot(normal_direction, v_avg) * (p_rr - p_ll)
end

@inline function flux_nonconservative_surface_artiano_td(p_ll, p_rr, v_interface,
                                                         equations::CompressibleEulerAtmo,
                                                         td_equation::EnergyInternal)
    return -v_interface * (p_rr - p_ll)
end
end # @muladd
