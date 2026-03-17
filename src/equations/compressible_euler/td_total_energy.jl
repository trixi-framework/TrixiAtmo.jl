@muladd begin
#! format: noindent

@doc raw"""
    TotalEnergy
"""
struct TotalEnergy{ThermodynamicStateType} <: AbstractThermodynamicEquation
    td_state::ThermodynamicStateType
end

varname_td(::typeof(cons2cons), ::TotalEnergy) = "rho_E"
varname_td(::typeof(cons2prim), ::TotalEnergy) = "p"

# TODO: T_ref!
@inline function temperature(u, equations::CompressibleEulerAtmo{NDIMS}, ::TotalEnergy,
                             td_state::IdealGas{RealType}) where
                            {NDIMS, RealType}
    @unpack cv = td_state

    rho = u[NDIMS+2]
    rho_moment = vars_moment(u, equations)

    # Kinetic energy times rho
    rho_Ekin = 0.5f0 * dot(rho_moment, rho_moment) / rho

    return (u[NDIMS+1] - rho_Ekin) / (rho * cv)
end

# TODO: T_ref!
@inline function temperature(u, equations::CompressibleEulerAtmo{NDIMS}, ::TotalEnergy,
                             td_state::IdealGasesAndLiquids{RealType}) where
                            {NDIMS, RealType}
    @unpack cv_gas, c_condens, gas_constant, latent_heat = td_state

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
    rho_cv = dot(rho_gas, cv_gas) + dot(rho_condens, c_condens)

    return (u[NDIMS+1] - rho_Ekin - rho_Elat) / rho_cv
end

@inline function pressure(u, equations::CompressibleEulerAtmo, ::TotalEnergy,
                          td_state::Union{IdealGas,IdealGasesAndLiquids})
    @unpack gas_constant = td_state
   
    rho_gas = vars_gas(u, equations)

    T = temperature(u, equations)

    # Gas constant of mixture times rho
    rho_R = sum(rho_gas .* gas_constant)

    return  rho_R * T
end

@inline function speed_of_sound(u, equations::CompressibleEulerAtmo, ::TotalEnergy,
                                td_state::IdealGasesAndLiquids)
    @unpack gas_constant = td_state
   
    rho_gas = vars_gas(u, equations)
    rho_condens = vars_condens(u, equations)

    cv_m = cv_total(rho_gas, rho_condens, td_state)
    R_m = R_total(rho_gas, td_state)
    cp_m = cv_m + R_m
    T = temperature(u, equations)
#@info T
    return sqrt(cp_m / cv_m * R_m * T)
end

@inline function flux_td(u, equations::CompressibleEulerAtmo{NDIMS}, ::TotalEnergy,
                         ::Union{IdealGas,IdealGasesAndLiquids}) where {NDIMS}

    p = pressure(u, equations)
    return u[NDIMS+1] + p
end

@inline function flux_lmars_td(p, v, ::CompressibleEulerAtmo{NDIMS, NVARS}, ::TotalEnergy) where {NDIMS, NVARS}
    return SVector(ntuple(i->0, Val(NDIMS))...,
                   p * v,
                   ntuple(i->0, Val(NVARS-NDIMS-1))...)
end

function prim2cons_td(prim, equations::CompressibleEulerAtmo{NDIMS}, ::TotalEnergy,
                      td_state::IdealGas) where {NDIMS}
    prim_moment = vars_moment(prim, equations)

    T = prim[NDIMS+1] / (prim[NDIMS+2] * td_state.gas_constant)
    # internal + kinetic energy
    return prim[NDIMS+2] * (td_state.cv * T +
                           0.5f0 * dot(prim_moment, prim_moment))
end

# TODO check
# TODO T_ref
function prim2cons_td(prim, equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS},
                      ::TotalEnergy, td_state::IdealGasesAndLiquids) where
                      {NDIMS, NVARS, NGAS}

    @unpack cv_gas, c_condens, gas_constant, latent_heat = td_state

    rho_total = prim2density_total(prim, equations)
    prim_gas = vars_gas(prim, equations)
    prim_vapor = var_vapor(prim, equations)
    prim_condens = vars_condens(prim, equations)
    prim_moment = vars_moment(prim, equations)

    # calculate temperature
    rho_Rm = gas_constant[1] * prim_gas[1] +
             dot(gas_constant[SVector{NGAS-1}(2:NGAS)],
                 prim_gas[SVector{NGAS-1}(2:NGAS)]) * rho_total
    T = prim[NDIMS+1] / rho_Rm

    # calculate energy
    Eint = (cv_gas[1] * prim_gas[1] / rho_total + 
            dot(cv_gas[SVector{NGAS-1}(2:NGAS)],
                prim_gas[SVector{NGAS-1}(2:NGAS)]) +
            dot(c_condens, prim_condens)) * T

    # latent heats
    Elatent = latent_heat * prim_vapor

    # kinetic energy
    Ekinetic = 0.5f0 * dot(prim_moment, prim_moment)

    return (Eint + Elatent + Ekinetic) * rho_total
end

end # @muladd