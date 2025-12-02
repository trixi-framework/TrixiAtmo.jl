@muladd begin
#! format: noindent

#
# These equations need to be compatible to Trixi.jl's interface
# therefore delegate anything special to thermodynamics etc.
# -> we need those as structs, therefore CompressibleEulerAtmo cannot be an abstract type
#

@doc raw"""
    CompressibleEulerAtmo

    type parameters:
    - NDIMS    dimension (2 or 3)
    - NVARS    number of total variables
    - NPASSIVE number of passive variables
               - not included in overall density
               - counted as part of NVARS

    members:
    - parameters    physical constants etc
    - td_state      thermodynamic state: used to ?
                      TODO: naming? needed?
    - td_equation   implementations for the "last conserved quantity" (energy, pot temp,..) 
                     TODO: naming?
    - microphysics  physical models for moisture

    conservative variables:
    - rho {u,v,w}  Cartesian momentum density components
    - rho TD       conserved thermodynamic quantitiy density
                   (given by thermodynamic_equation)
    - rho_X        density of all species (active and passive)

    rho (total density) = sum_A rho_X (all active species) 

    primitive variables:
    - {u,v,w}      Cartesian velocities components
    - TD           conserved thermodynamic quantitiy
                   (given by thermodynamic_equation)
    - rho_d        density of first species (e.g. dry air)
    - r_X          fraction w.r.t. total density of further species
"""
struct CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE,
                          ParametersType,
                          ThermodynamicStateType,
                          ThermodynamicEquationType,
                          MicrophysicsType} <: AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}

    parameters::ParametersType
    td_state::ThermodynamicStateType
    td_equation::ThermodynamicEquationType
    microphysics::MicrophysicsType

    function CompressibleEulerAtmo{NDIMS}(; parameters,
        thermodynamic_state::AbstractThermodynamicState{RealType, NGAS, NCONDENS, NPRECIP},
        thermodynamic_equation,
        microphysics,
        n_vars_passive = 0) where {NDIMS, RealType, NGAS, NCONDENS, NPRECIP}

        # gaseous components + airborn liquid components + precipitating components
        # + passive components + momentum equations (NDIMS) + thermodynamic equation
        n_vars = NGAS + NCONDENS + NPRECIP + n_vars_passive + NDIMS + 1

        return new{NDIMS, n_vars, NGAS, NCONDENS, NPRECIP, n_vars_passive,
                   typeof(parameters), typeof(thermodynamic_state),
                   typeof(thermodynamic_equation), typeof(microphysics)}(
                   parameters, thermodynamic_state, thermodynamic_equation, microphysics)
    end
end

have_nonconservative_terms(equations::CompressibleEulerAtmo) =
    have_nonconservative_terms(equations.td_equation)

@inline function Base.real(::CompressibleEulerAtmo)
    return Base.real(equations.td_state)
end

@inline function vars_moment(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[SVector{NDIMS}(1:NDIMS)]
end

@inline function vars_gas(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[SVector{NGAS}(NDIMS+2:NDIMS+1+NGAS)]
end

@inline function var_vapor(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, 1, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NCONDENS, NPRECIP, NPASSIVE}
    return 0
end

@inline function var_vapor(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[NDIMS+3]
end

@inline function vars_condens(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[SVector{NCONDENS}(NDIMS+2+NGAS:NDIMS+1+NGAS+NCONDENS)]
end

@inline function vars_airborn(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[SVector{NGAS+NCONDENS}(NDIMS+2:NDIMS+1+NGAS+NCONDENS)]
end

@inline function vars_precip(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[SVector{NCONDENS}(NDIMS+2+NGAS+NCONDENS:NDIMS+1+NGAS+NCONDENS+NPRECIP)]
end

@inline function vars_passive(u,
    ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return u[SVector{NPASSIVE}(NVARS-NPASSIVE+1:end)]
end

function varnames(variables::typeof(cons2cons),
    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return (ntuple(i -> "rho_v$i", Val(NDIMS))...,
            varname_td(variables, equations.td_equation),
            ("rho" .* varnames_gas(variables, equations.td_state))...,
            ("rho" .* varnames_liquid(variables, equations.td_state))...,
            ("rho" .* varnames_precip(variables, equations.td_state))...,
            ntuple(i -> "rho_chi_$i", Val(NPASSIVE))...)
end

function varnames(variables::typeof(cons2prim),
    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}) where 
    {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    return (ntuple(i -> "v$i", Val(NDIMS))...,
            varname_td(variables, equations.td_equation),
            varnames_gas(variables, equations.td_state)...,
            varnames_liquid(variables, equations.td_state)...,
            varnames_precip(variables, equations.td_state)...,
            ntuple(i -> "chi_$i", Val(NPASSIVE))...)
end

@inline function velocities_precip(u, equations::CompressibleEulerAtmo)
    @error "velocities_precip not implemented for more than one species"
end

@inline function velocities_precip(u,
    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, 0}) where {NDIMS, NVARS, NGAS}
    return SVector{0, eltype(u)}()
end

# TODO: currently rain only
@inline function velocities_precip(u,
    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, 1}) where {NDIMS, NVARS, NGAS}

    return velocity_rain(u, equations, equations.microphysics)
end

# Convert full vector of primitive variables to conservative quantities
@inline function prim2cons(prim,
                           equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {NDIMS, NVARS}
    rho = prim2density_total(prim, equations)
    return(SVector(ntuple(i -> rho * prim[i], Val(NDIMS))...,
                   prim2cons_td(prim, equations, equations.td_equation, equations.td_state),
                   prim[NDIMS+2],
                   ntuple(i -> rho * prim[NDIMS+2+i], Val(NVARS-NDIMS-2))...))
end

# Convert full vector of conservative variables to primitive quantities
@inline function cons2prim(cons,
                           equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {NDIMS, NVARS}
    rho = density_total(cons, equations)
    return(SVector(ntuple(i -> cons[i] / rho, Val(NDIMS))...,
                   pressure(cons, equations),
                   cons[NDIMS+2],
                   ntuple(i -> cons[NDIMS+2+i] / rho, Val(NVARS-NDIMS-2))...))
end

# Temperature needs to be calculated depending on the thermodynamic equation
@inline function temperature(cons, equations::CompressibleEulerAtmo)
    return temperature(cons, equations, equations.td_equation, equations.td_state)
end

# Pressure needs to be calculated depending on the thermodynamic equation
@inline function pressure(cons, equations::CompressibleEulerAtmo)
    return pressure(cons, equations, equations.td_equation, equations.td_state)
end

# Pressure needs to be calculated depending on the thermodynamic equation
@inline function speed_of_sound(cons, equations::CompressibleEulerAtmo)
    return speed_of_sound(cons, equations, equations.td_equation, equations.td_state)
end

# Convert conservative variables to entropy
# TODO
@inline function cons2entropy(cons, equations::CompressibleEulerAtmo)

    rho_total = density_total(cons, equations)
    rho_gas = vars_gas(cons, equations)
    rho_liquid = vars_gas(cons, equations)
    T = 270.0 # temperature(cons, equations)
    
    #s_gas = entropies_gas(rho_gas, T, equations.td_state)
    #s_liquid = entropies_liquid(rho_liquid, T, equations.td_state)
    
  return cons
end

# Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer,
                      equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    normal_direction = SVector{NDIMS}(i == orientation ? 1 : 0 for i in 1:NDIMS)
    return flux(u, normal_direction, equations)
end

@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS}) where
                      {NDIMS, NVARS, NGAS, NCONDENS}
    rho_total = density_total(u, equations)
    v_normal = dot(vars_moment(u, equations), normal_direction) / rho_total
    p = pressure(u, equations)

    # momentum equations
    f_mom = vars_moment(u, equations) .* v_normal + normal_direction .* p

    # thermodynamic_equation
    f_td = flux_td(u, equations, equations.td_equation, equations.td_state) * v_normal

    # mass equations
    f_mass_air = vars_airborn(u, equations) .* v_normal
    f_mass_precip1 = vars_precip(u, equations) .* v_normal
    f_mass_passive = vars_passive(u, equations) .* v_normal

    f_mom_precip, f_td_precip, f_mass_precip2 = flux_precip(u, normal_direction, equations)

    ret = SVector((f_mom + f_mom_precip)..., f_td + f_td_precip,
                    f_mass_air..., (f_mass_precip1 + f_mass_precip2)..., f_mass_passive...)
    
    return ret
end

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {NDIMS, NVARS}
    a = flux_lmars.speed_of_sound
    norm_ = norm(normal_direction)

    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_n_ll = dot(vars_moment(u_ll, equations), normal_direction) / rho_total_ll
    v_n_rr = dot(vars_moment(u_rr, equations), normal_direction) / rho_total_rr
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    rho = 0.5f0 * (rho_total_ll + rho_total_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_n_rr - v_n_ll) / norm_
    v_interface = 0.5f0 * (v_n_ll + v_n_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f = u_ll * v_interface
    else
        f = u_rr * v_interface
    end

    # additional terms in momentum equation, pad with zeros
    f_mom = SVector((normal_direction * p_interface)...,
                    ntuple(i->0, Val(NVARS-NDIMS))...)

    return f + f_mom
end

@inline function flux_precip(u, normal_direction::AbstractVector,
    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS}) where
                      {NDIMS, NVARS, NGAS, NCONDENS}

    # precipitation along third direction
    # TODO this is not correct for curved geometries (earth)
    # TODO earth normal componente instead of normal_direction[NDIMS]

    rho_total = density_total(u, equations)
    v_precip = velocities_precip(u, equations)
    f_mom_z = vars_moment(u, equations)[NDIMS] * normal_direction[NDIMS] / rho_total *
              dot(v_precip, vars_precip(u, equations))
    f_mass_precip = vars_precip(u, equations) .* v_precip * normal_direction[NDIMS]

    return SVector(ntuple(i -> 0, Val(NDIMS-1))..., f_mom_z),
           0, 
           f_mass_precip
end


# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)

    v_ll = dot(vars_moment(u_ll, equations), normal_direction) / rho_total_ll
    v_rr = dot(vars_moment(u_rr, equations), normal_direction) / rho_total_rr

    c_ll = speed_of_sound(u_ll, equations)
    c_rr = speed_of_sound(u_rr, equations)

    # TODO: assumes z axis normal to surface
    v_precip_ll = sum(abs.(velocities_precip(u_ll, equations) .*
        normal_direction[NDIMS]))
    v_precip_rr = sum(abs.(velocities_precip(u_rr, equations) .*
        normal_direction[NDIMS]))

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + abs(v_precip_ll) + c_ll * norm_,
               abs(v_rr) + abs(v_precip_rr) + c_rr * norm_)
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    v = velocity(u, equations)
    c = speed_of_sound(u, equations)                           

    # TODO: assumed z direction
    v_precip_z = maximum(abs, velocities_precip(u, equations))
    v_precip = SVector{NDIMS}(i == NDIMS ? v_precip_z : 0 for i in 1:NDIMS)

    return abs.(v) .+ c .+ v_precip
end

"""
    boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,
                                 equations::CompressibleEulerAtmo)

Determine the boundary numerical surface flux for a slip wall condition.
Imposes a zero normal velocity at the wall.
Density is taken from the internal solution state and pressure is computed as an
exact solution of a 1D Riemann problem. Further details about this boundary state
are available in the paper:
- J. J. W. van der Vegt and H. van der Ven (2002)
  Slip flow boundary conditions in discontinuous Galerkin discretizations of
  the Euler equations of gas dynamics
  [PDF](https://reports.nlr.nl/bitstream/handle/10921/692/TP-2002-300.pdf?sequence=1)

Details about the 1D pressure Riemann solution can be found in Section 6.3.3 of the book
- Eleuterio F. Toro (2009)
  Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
  3rd edition
  [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)

Should be used together with [`UnstructuredMesh2D`](@ref), [`P4estMesh`](@ref), or [`T8codeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
    x, t, surface_flux_function,
    equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {NDIMS, NVARS}

    rho_gas = vars_gas(u_inner, equations)
    rho_condens = vars_condens(u_inner, equations)

    # total quantities
    rho = density_total(u_inner, equations)
    gamma = gamma_total(rho_gas, rho_condens, equations.td_state)

    # normal velocity
    normal = normal_direction / norm(normal_direction)
    v_normal = dot(vars_moment(u_inner, equations), normal) / rho
    
    # compute pressure
    p_local = pressure(u_inner, equations)

    # Get the solution of the pressure Riemann problem
    if v_normal <= 0
        sound_speed = sqrt(gamma * p_local / rho) # local sound speed
        p_star = p_local *
                 (1 + 0.5f0 * (gamma - 1) * v_normal / sound_speed)^(2 *
                                                                     gamma /
                                                                     (gamma - 1.0f0))
    else # v_normal > 0
        A = 2 / ((gamma + 1) * rho)
        B = p_local * (gamma - 1) / (gamma + 1)
        p_star = p_local +
                 0.5f0 * v_normal / A *
                 (v_normal + sqrt(v_normal^2 + 4 * A * (p_local + B)))
    end

    # For the slip wall we directly set the flux as the normal velocity is zero
    return SVector((normal_direction .* p_star)...,
                   ntuple(i->0, Val(NVARS-NDIMS))...)
end

"""
    boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                 surface_flux_function, equations::CompressibleEulerAtmo)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
    surface_flux_function, equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}

    normal_direction = SVector{NDIMS}(i == orientation ? 1 : 0 for i in 1:NDIMS)

    # compute and return the flux using `boundary_condition_slip_wall` routine above
    return boundary_condition_slip_wall(u_inner, normal_direction, direction,
                                        x, t, surface_flux_function, equations)
end

"""
    boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t,
                                 surface_flux_function, equations::CompressibleEulerEquations3D)

Should be used together with [`StructuredMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                              direction, x, t,
                                              surface_flux_function,
                                              equations::CompressibleEulerAtmo)
    # flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
    # to be inward pointing on the -x, -y, and -z sides due to the orientation convention used by StructuredMesh
    if isodd(direction)
        boundary_flux = -boundary_condition_slip_wall(u_inner, -normal_direction,
                                                      x, t, surface_flux_function,
                                                      equations)
    else
        boundary_flux = boundary_condition_slip_wall(u_inner, normal_direction,
                                                     x, t, surface_flux_function,
                                                     equations)
    end

    return boundary_flux
end

end # @muladd
