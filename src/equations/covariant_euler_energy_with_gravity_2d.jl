@muladd begin
#! format: noindent

@doc raw"""
    CovariantEulerEnergyEquationsWithGravity2D{GlobalCoordinateSystem} <:  
        AbstractCovariantEquations{2, 2, GlobalCoordinateSystem, 4}
## References
- Comparison of two Euler equation sets in a Discontinuous Galerkin
solver for atmospheric modelling (BRIDGE v0.9)
Michael Baldauf and Florian Prill
"""
struct CovariantEulerEnergyEquationsWithGravity2D{GlobalCoordinateSystem,
                                                  RealT <: Real} <:
       AbstractCovariantEulerEquations{2, GlobalCoordinateSystem, 4}
    gamma::RealT  # 
    inv_gamma_minus_one::RealT  # 
    global_coordinate_system::GlobalCoordinateSystem
    function CovariantEulerEnergyEquationsWithGravity2D(gamma::RealT;
                                                        global_coordinate_system = GlobalCartesianCoordinates()) where {RealT <:
                                                                                                                        Real}
        return new{typeof(global_coordinate_system), RealT}(gamma, inv(gamma - 1))
    end
end

have_nonconservative_terms(::CovariantEulerEnergyEquationsWithGravity2D) = False()

# The conservative variables are the height and contravariant momentum components
function varnames(::typeof(cons2cons), ::CovariantEulerEnergyEquationsWithGravity2D)
    return ("rho", "rho_vcon1", "rho_vcon2", "rho_e")
end

# The primitive variables are the height and contravariant velocity components
function varnames(::typeof(cons2prim), ::CovariantEulerEnergyEquationsWithGravity2D)
    return ("rho", "vcon1", "vcon2", "p")
end

# The change of variables contravariant2global converts the two local contravariant vector 
# components u[2] and u[3] to the three global vector components specified by 
# equations.global_coordinate_system (e.g. spherical or Cartesian). This transformation 
# works for both primitive and conservative variables, although varnames refers 
# specifically to transformations from conservative variables.
function varnames(::typeof(contravariant2global),
                  ::CovariantEulerEnergyEquationsWithGravity2D)
    return ("rho", "v1", "v2", "e")
end

# Convenience functions to extract physical variables from state vector
@inline density(u, ::CovariantEulerEnergyEquationsWithGravity2D) = u[1]

@inline velocity_contravariant(u, ::CovariantEulerEnergyEquationsWithGravity2D) = SVector(u[2] /
                                                                                          u[1],
                                                                                          u[3] /
                                                                                          u[1])
@inline momentum_contravariant(u, ::CovariantEulerEnergyEquationsWithGravity2D) = SVector(u[2],
                                                                                          u[3])

@inline total_energy(u, ::CovariantEulerEnergyEquationsWithGravity2D) = u[4]

@inline energy_density(u, ::CovariantEulerEnergyEquationsWithGravity2D) = u[4] / u[1]

@inline function kinetic_energy(u, aux_vars, equations::AbstractCovariantEulerEquations)
    rho = density(u, equations)
    vcon = velocity_contravariant(u, equations)
    Gcov = metric_covariant(aux_vars, equations)
    return 0.5f0 * dot(Gcov * vcon, vcon) * rho
end

@inline function pressure(u, aux_vars, equations::AbstractCovariantEulerEquations)
    rho = density(u, equations)
    E_total = total_energy(u, equations)
    ekin = kinetic_energy(u, aux_vars, equations)
    phi = geopotential(aux_vars, equations)
    return (equations.gamma - 1) * (E_total - ekin - rho * phi)
end

@inline function cons2prim(u, aux_vars,
                           equations::CovariantEulerEnergyEquationsWithGravity2D)
    rho = density(u, equations)
    vcon = velocity_contravariant(u, equations)
    p = pressure(u, aux_vars, equations)
    return SVector(rho, vcon[1], vcon[2], p)
end

@inline function prim2cons(u, aux_vars,
                           equations::CovariantEulerEnergyEquationsWithGravity2D)
    rho, vcon1, vcon2, p = u
    vcon = SVector(vcon1, vcon2)
    Gcov = metric_covariant(aux_vars, equations)
    ekin = 0.5f0 * dot(Gcov * vcon, vcon) * rho
    phi = geopotential(aux_vars, equations)
    rho_e_total = p * equations.inv_gamma_minus_one + ekin + rho * phi
    return SVector(rho, rho * vcon1, rho * vcon2, rho_e_total)
end

@inline function cons2entropy(u, aux_vars,
                              equations::CovariantEulerEnergyEquationsWithGravity2D)
    Gcov = metric_covariant(aux_vars, equations)
    rho, rho_vcon1, rho_vcon2, rho_e_total = u
    p = pressure(u, aux_vars, equations)
    vcon = SVector(rho_vcon1 / rho, rho_vcon2 / rho)
    vcov = Gcov * vcon
    s = log(p / rho^equations.gamma)

    # Entropy variables (covariant form)
    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         0.5f0 * rho / p * dot(vcov, vcon)
    w2 = rho / p * vcov[1]
    w3 = rho / p * vcov[2]
    w4 = -rho / p

    return SVector(w1, w2, w3, w4)
end

# Convert contravariant momentum components to the global coordinate system
@inline function contravariant2global(u, aux_vars,
                                      equations::CovariantEulerEnergyEquationsWithGravity2D)
    v1, v2 = basis_covariant(aux_vars, equations) *
             SVector(u[2], u[3])
    return SVector(u[1], v1, v2, u[4])
end

# Convert momentum components in the global coordinate system to contravariant components
@inline function global2contravariant(u, aux_vars,
                                      equations::CovariantEulerEnergyEquationsWithGravity2D)
    vcon1, vcon2 = basis_contravariant(aux_vars, equations) *
                   SVector(u[2], u[3])
    return SVector(u[1], vcon1, vcon2, u[4])
end

# Flux as a function of the state vector u, as well as the auxiliary variables aux_vars, 
# which contain the geometric information required for the covariant form
@inline function flux(u, aux_vars, orientation::Integer,
                      equations::AbstractCovariantEulerEquations)
    # Geometric variables
    Gcon = metric_contravariant(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    rho_vcon = momentum_contravariant(u, equations)
    vcon = velocity_contravariant(u, equations)
    E_total = total_energy(u, equations)

    # Compute and store the pressure and Energy momentum tensor components in the desired orientation
    p = pressure(u, aux_vars, equations)
    T = rho_vcon * vcon[orientation] + p * Gcon[:, orientation]

    return SVector(J * rho_vcon[orientation], J * T...,
                   J * vcon[orientation] * (E_total + p))
end

@inline function flux(u, aux_vars, normal_direction::AbstractVector,
                      equations::AbstractCovariantEulerEquations)
    # Geometric variables
    Gcon = metric_contravariant(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    rho = density(u, equations)
    rho_vcon = momentum_contravariant(u, equations)
    E_total = total_energy(u, equations)

    # Compute and store the pressure and Energy momentum tensor components in the desired direction
    vcon = dot(rho_vcon, normal_direction) / rho
    p = pressure(u, aux_vars, equations)
    T = rho_vcon * vcon + p * (Gcon * normal_direction)

    return SVector(J * dot(rho_vcon, normal_direction), J * T...,
                   J * vcon * (E_total + p))
end

# Maximum wave speeds with respect to the covariant basis
@inline function max_abs_speeds(u, aux_vars,
                                equations::CovariantEulerEnergyEquationsWithGravity2D)
    rho, vcon1, vcon2, p = cons2prim(u, aux_vars, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    c1 = sqrt(equations.gamma * p / rho * Gcon[1, 1])
    c2 = sqrt(equations.gamma * p / rho * Gcon[2, 2])

    return abs(vcon1) + c1, abs(vcon2) + c2
end

# Maximum wave speed along the normal direction in reference space
@inline function max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               normal_direction::AbstractVector,
                               equations::CovariantEulerEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, aux_vars_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, aux_vars_rr, equations)

    # Calcualte the normal velocities and sound speeds
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll *
                (normal_direction' * metric_contravariant(aux_vars_ll, equations) *
                 normal_direction))

    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr *
                (normal_direction' * metric_contravariant(aux_vars_rr, equations) *
                 normal_direction))

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

@inline function boundary_condition_slip_wall(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function,
                                              equations::CovariantEulerEnergyEquationsWithGravity2D)
    # Reflect the normal contravariant velocity component and keep the tangential component unchanged
    u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4])

    # Compute the flux at the boundary using the surface flux function
    flux = surface_flux_function(u_inner, u_boundary, aux_inner, aux_inner,
                                 normal_direction, equations)
    return flux
end

function source_terms_gravity(u, x, t, aux_vars,
                              equations::CovariantEulerEnergyEquationsWithGravity2D)
    rho = u[1]
    Gcon = metric_contravariant(aux_vars, equations)
    Gcov = metric_covariant(aux_vars, equations)
    return SVector(0.0, 0.0, -9.81 * rho * Gcon[2, 2] * sqrt(Gcov[2, 2]), 0.0)
end
end # @muladd
