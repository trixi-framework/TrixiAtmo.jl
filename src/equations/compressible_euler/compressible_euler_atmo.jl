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
    - NGAS     number of gaseous species
    - NCONDENS number of condensed species
    - NPRECIP  number of precipitating species
    - NPASSIVE number of passive variables
               - not included in overall density
               - counted as part of NVARS
    - NAUX     number of additional auxiliary variables

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
struct CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX,
                             ParametersType,
                             ThermodynamicStateType,
                             ThermodynamicEquationType,
                             MicrophysicsType} <:
       AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}
    parameters::ParametersType
    td_state::ThermodynamicStateType
    td_equation::ThermodynamicEquationType
    microphysics::MicrophysicsType

    function CompressibleEulerAtmo(; n_dims, parameters,
                                   thermodynamic_state::AbstractThermodynamicState{RealType,
                                                                                   NGAS,
                                                                                   NCONDENS,
                                                                                   NPRECIP},
                                   thermodynamic_equation,
                                   microphysics = nothing,
                                   n_vars_passive = 0,
                                   n_vars_aux = 0,
                                   NAUX = false) where {RealType, NGAS, NCONDENS,
                                                        NPRECIP}
        n_vars = NGAS + NCONDENS + NPRECIP + n_vars_passive + n_dims + 1

        return new{n_dims, n_vars, NGAS, NCONDENS, NPRECIP, n_vars_passive, n_vars_aux,
                   typeof(parameters), typeof(thermodynamic_state),
                   typeof(thermodynamic_equation), typeof(microphysics)}(parameters,
                                                                         thermodynamic_state,
                                                                         thermodynamic_equation,
                                                                         microphysics)
    end
end

# TODO: at the moment aux vars means gravity means noncons terms
@inline have_nonconservative_terms(::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}) where {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX} = static(NAUX >
                                                                                                                                                                                          0)

@inline Trixi.have_aux_node_vars(::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}) where {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX} = static(NAUX >
                                                                                                                                                                                        0)

@inline Trixi.n_aux_node_vars(::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}) where {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX} = NAUX

@inline function Base.real(::CompressibleEulerAtmo)
    return Base.real(equations.td_state)
end

@inline function vars_moment(u,
                             ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                     NPRECIP, NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[SVector{NDIMS}(1:NDIMS)]
end

@inline function var_td(u,
                        ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP,
                                                NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[NDIMS + 1]
end

@inline function vars_gas(u,
                          ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS, NPRECIP,
                                                  NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[SVector{NGAS}((NDIMS + 2):(NDIMS + 1 + NGAS))]
end

@inline function var_vapor(u,
                           ::CompressibleEulerAtmo{NDIMS, NVARS, 1, NCONDENS, NPRECIP,
                                                   NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return 0
end

@inline function var_vapor(u,
                           ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                   NPRECIP, NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[NDIMS + 3]
end

@inline function vars_condens(u,
                              ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                      NPRECIP, NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[SVector{NCONDENS}((NDIMS + 2 + NGAS):(NDIMS + 1 + NGAS + NCONDENS))]
end

@inline function vars_airborn(u,
                              ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                      NPRECIP, NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[SVector{NGAS + NCONDENS}((NDIMS + 2):(NDIMS + 1 + NGAS + NCONDENS))]
end

@inline function vars_precip(u,
                             ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                     NPRECIP, NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[SVector{NPRECIP}((NDIMS + 2 + NGAS + NCONDENS):(NDIMS + 1 + NGAS + NCONDENS + NPRECIP))]
end

@inline function vars_passive(u,
                              ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                      NPRECIP, NPASSIVE, NAUX}) where
                 {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return u[SVector{NPASSIVE}((NVARS - NPASSIVE + 1):end)]
end

function varnames(variables::typeof(cons2cons),
                  equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                   NPRECIP, NPASSIVE, NAUX}) where
         {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return (ntuple(i -> "rho_v$i", Val(NDIMS))...,
            varname_td(variables, equations.td_equation),
            ("rho_" .* varnames_gas(variables, equations.td_state))...,
            ("rho_" .* varnames_liquid(variables, equations.td_state))...,
            ("rho_" .* varnames_precip(variables, equations.td_state))...,
            ntuple(i -> "rho_chi_$i", Val(NPASSIVE))...)
end

function varnames(variables::typeof(cons2prim),
                  equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                                   NPRECIP, NPASSIVE, NAUX}) where
         {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE, NAUX}
    return (ntuple(i -> "v$i", Val(NDIMS))...,
            varname_td(variables, equations.td_equation),
            (SVector{NGAS}("rho_", ntuple(i -> "r_", Val(NGAS - 1))...) .*
             varnames_gas(variables, equations.td_state))...,
            ("r_" .* varnames_liquid(variables, equations.td_state))...,
            ("r_" .* varnames_precip(variables, equations.td_state))...,
            ntuple(i -> "chi_$i", Val(NPASSIVE))...)
end

function varnames(::typeof(Trixi.cons2aux),
                  ::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS,
                                          NPRECIP, NPASSIVE, 1}) where
         {NDIMS, NVARS, NGAS, NCONDENS, NPRECIP, NPASSIVE}
    ("geopotential",)
end

@inline function velocities_precip(u, equations::CompressibleEulerAtmo)
    @error "velocities_precip not implemented for more than one species"
end

@inline function velocities_precip(u,
                                   equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS,
                                                                    0}) where {NDIMS,
                                                                               NVARS,
                                                                               NGAS}
    return SVector{0, eltype(u)}()
end

# TODO: currently rain only
@inline function velocities_precip(u,
                                   equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS,
                                                                    1}) where {NDIMS,
                                                                               NVARS,
                                                                               NGAS}
    return velocity_rain(u, equations, equations.microphysics)
end

# Convert full vector of primitive variables to conservative quantities
@inline function prim2cons(prim,
                           equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {NDIMS,
                                                                                  NVARS}
    rho = prim2density_total(prim, equations)
    return (SVector(ntuple(i -> rho * prim[i], Val(NDIMS))...,
                    prim2cons_td(prim, equations, equations.td_equation,
                                 equations.td_state),
                    prim[NDIMS + 2],
                    ntuple(i -> rho * prim[NDIMS + 2 + i], Val(NVARS - NDIMS - 2))...))
end

# Convert full vector of conservative variables to primitive quantities
@inline cons2prim(cons, aux, equations::CompressibleEulerAtmo) = cons2prim(cons,
                                                                           equations)
@inline function cons2prim(cons,
                           equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {NDIMS,
                                                                                  NVARS}
    rho = density_total(cons, equations)
    return (SVector(ntuple(i -> cons[i] / rho, Val(NDIMS))...,
                    pressure(cons, equations),
                    cons[NDIMS + 2],
                    ntuple(i -> cons[NDIMS + 2 + i] / rho, Val(NVARS - NDIMS - 2))...))
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
@inline cons2entropy(cons, aux, equations::CompressibleEulerAtmo) = cons2entropy(cons,
                                                                                 equations)
@inline function cons2entropy(cons, equations::CompressibleEulerAtmo)
    return cons2entropy(cons, equations, equations.td_equation, equations.td_state)
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
    #    v_precip_ll = sum(abs.(velocities_precip(u_ll, equations) .*
    #        normal_direction[NDIMS]))
    #   v_precip_rr = sum(abs.(velocities_precip(u_rr, equations) .*
    #       normal_direction[NDIMS]))

    norm_ = norm(normal_direction)
    #+ abs(v_precip_ll)
    #+ abs(v_precip_rr)
    return max(abs(v_ll) + c_ll * norm_,
               abs(v_rr) + c_rr * norm_)
end

@inline max_abs_speeds(u, aux, equations::CompressibleEulerAtmo) = max_abs_speeds(u,
                                                                                  equations)
@inline function max_abs_speeds(u,
                                equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    v = velocity(u, equations)
    c = speed_of_sound(u, equations)

    # TODO: assumed z direction
    #v_precip_z = maximum(abs, velocities_precip(u, equations))
    #v_precip = SVector{NDIMS}(i == NDIMS ? v_precip_z : 0 for i in 1:NDIMS)

    return abs.(v) .+ c #.+ v_precip
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
                                              equations::CompressibleEulerAtmo{NDIMS,
                                                                               NVARS}) where {
                                                                                              NDIMS,
                                                                                              NVARS
                                                                                              }
    rho_gas = vars_gas(u_inner, equations)
    rho_condens = vars_condens(u_inner, equations)

    # total quantities
    rho = density_total(u_inner, equations)
    gamma = gamma_total(rho_gas, rho_condens, equations.td_state)

    # normal velocity
    v_normal = dot(vars_moment(u_inner, equations), normal_direction) /
               (rho * norm(normal_direction))

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
    return SVector((normal_direction * p_star)...,
                   ntuple(i -> 0, Val(NVARS - NDIMS))...)
end

"""
    boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                 surface_flux_function, equations::CompressibleEulerAtmo)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                              surface_flux_function,
                                              equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
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

@inline function boundary_condition_slip_wall_simple(u_inner, orientation, direction, x,
                                                     t, surface_flux_function,
                                                     equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    # flip the normal component of momentum densities
    u_mom = vars_moment(u_inner, equations)
    u_boundary_mom = SVector{NDIMS}(i == orientation ? -u_mom[i] : u_mom[i]
                                    for i in 1:NDIMS)
    u_boundary = SVector(u_boundary_mom...,
                         var_td(u_inner, equations),
                         vars_airborn(u_inner, equations)...,
                         vars_precip(u_inner, equations)...,
                         vars_passive(u_inner, equations)...)

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

@inline function boundary_condition_slip_wall_simple(u_inner, aux_inner,
                                                     normal_direction::AbstractVector,
                                                     x, t,
                                                     surface_flux_functions::Tuple,
                                                     equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # compute the normal velocity
    u_mom = vars_moment(u_inner, equations)
    u_mom_normal = dot(vars_moment(u_inner, equations), normal)
    u_boundary_mom = u_mom - 2 * u_mom_normal * normal

    # create the "external" boundary solution state
    u_boundary = SVector(u_boundary_mom...,
                         var_td(u_inner, equations),
                         vars_airborn(u_inner, equations)...,
                         vars_precip(u_inner, equations)...,
                         vars_passive(u_inner, equations)...)

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, aux_inner, aux_inner,
                                 normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, aux_inner,
                                                 aux_inner, normal_direction, equations)
    return flux, noncons_flux
end

@inline function boundary_condition_slip_wall_simple(u_inner,
                                                     normal_direction::AbstractVector,
                                                     x, t,
                                                     surface_flux_functions::Tuple,
                                                     equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # compute the normal velocity
    u_mom = vars_moment(u_inner, equations)
    u_mom_normal = dot(vars_moment(u_inner, equations), normal)
    u_boundary_mom = u_mom - 2 * u_mom_normal * normal

    # create the "external" boundary solution state
    u_boundary = SVector(u_boundary_mom...,
                         var_td(u_inner, equations),
                         vars_airborn(u_inner, equations)...,
                         vars_precip(u_inner, equations)...,
                         vars_passive(u_inner, equations)...)

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

# Transform source terms, see abstract
# special case: no aux vars
function transform_source_terms_sum(source_terms::Tuple,
                                    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS,
                                                                     NCONDENS, NPRECIP,
                                                                     NPASSIVE, 0}) where {
                                                                                          NDIMS,
                                                                                          NVARS,
                                                                                          NGAS,
                                                                                          NCONDENS,
                                                                                          NPRECIP,
                                                                                          NPASSIVE
                                                                                          }
    function source_terms_transformed(u, x, t,
                                      equations::CompressibleEulerAtmo)
        # evaluate and add source terms
        source_ref = mapreduce(+, source_terms) do f
            ret = f(u, x, t, equations)
            SVector(ret..., ntuple(i -> 0, Val(NVARS - length(ret)))...)
        end
        return source_ref
    end
    return source_terms_transformed
end

# general case: with aux vars
function transform_source_terms_sum(source_terms::Tuple,
                                    equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS,
                                                                     NCONDENS, NPRECIP,
                                                                     NPASSIVE, NAUX}) where {
                                                                                             NDIMS,
                                                                                             NVARS,
                                                                                             NGAS,
                                                                                             NCONDENS,
                                                                                             NPRECIP,
                                                                                             NPASSIVE,
                                                                                             NAUX
                                                                                             }
    function source_terms_transformed(u, aux, x, t,
                                      equations::CompressibleEulerAtmo)
        # evaluate and add source terms
        source_ref = mapreduce(+, source_terms) do f
            ret = f(u, aux, x, t, equations)
            SVector(ret..., ntuple(i -> 0, Val(NVARS - length(ret)))...)
        end
        return source_ref
    end
    return source_terms_transformed
end
end # @muladd
