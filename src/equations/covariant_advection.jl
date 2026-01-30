@muladd begin
#! format: noindent

@doc raw"""
    CovariantLinearAdvectionEquation2D{GlobalCoordinateSystem} <:  
        AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}
Denoting the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) by 
$\nabla_j$ and summing over repeated indices, a variable-coefficient linear advection equation can be defined on a two-dimensional manifold in three-dimensional ambient space as
```math
\partial_t h + \nabla_j (hv^j) = 0.
```
We treat this problem as a system of equations in which the first variable is the scalar 
conserved quantity $h$, and the second two are the contravariant components $v^1$ and $v^2$ 
used in the expansion 
```math
\vec{v} = v^i \vec{a}_i =  v^1 \vec{a}_1 + v^2 \vec{a}_2,
```
where $\vec{a}_1 = \partial \vec{x} / \partial \xi^1$ and 
$\vec{a}_2 = \partial \vec{x} / \partial \xi^2$ are the so-called covariant basis vectors, 
and $\xi^1$ and $\xi^2$ are the local reference space coordinates. The velocity components 
are spatially varying but assumed to be constant in time, so we do not apply any flux or 
dissipation to such variables. The resulting system is then given on the reference element 
as 
```math
J \frac{\partial}{\partial t}
\left[\begin{array}{c} h \\ v^1 \\ v^2 \end{array}\right] 
+
\frac{\partial}{\partial \xi^1} 
\left[\begin{array}{c} J h v^1 \\ 0 \\ 0 \end{array}\right]
+ 
\frac{\partial}{\partial \xi^2} 
\left[\begin{array}{c} J h v^2 \\ 0 \\ 0 \end{array}\right] 
= 
\left[\begin{array}{c} 0 \\ 0 \\ 0 \end{array}\right],
```
where $J = \lVert\vec{a}^1 \times \vec{a}^2 \rVert$ is the area element. Note that the 
variable advection velocity components could alternatively be stored as auxiliary 
variables, similarly to the geometric information.
"""
struct CovariantLinearAdvectionEquation2D{GlobalCoordinateSystem} <:
       AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}
    global_coordinate_system::GlobalCoordinateSystem
    function CovariantLinearAdvectionEquation2D(;
                                                global_coordinate_system = GlobalCartesianCoordinates())
        return new{typeof(global_coordinate_system)}(global_coordinate_system)
    end
end

# The conservative variables are the scalar conserved quantity and two contravariant 
# velocity components.
function varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation2D)
    return ("h", "vcon1", "vcon2")
end

# Convenience function to extract the velocity
function velocity_contravariant(u, ::CovariantLinearAdvectionEquation2D)
    return SVector(u[2], u[3])
end

# Convert contravariant velocity components to the global coordinate system
@inline function contravariant2global(u, aux_vars,
                                      equations::CovariantLinearAdvectionEquation2D)
    v1, v2, v3 = basis_covariant(aux_vars, equations) *
                 velocity_contravariant(u, equations)
    return SVector(u[1], v1, v2, v3)
end

# Convert velocity components in the global coordinate system to contravariant components
@inline function global2contravariant(u, aux_vars,
                                      equations::CovariantLinearAdvectionEquation2D)
    vcon1, vcon2 = basis_contravariant(aux_vars, equations) * SVector(u[2], u[3], u[4])
    return SVector(u[1], vcon1, vcon2)
end

# Scalar conserved quantity and three global velocity components
function varnames(::typeof(contravariant2global),
                  ::CovariantLinearAdvectionEquation2D)
    return ("h", "v1", "v2", "v3")
end

# We will define the "entropy variables" here to just be the scalar variable in the first 
# slot, with zeros in all other positions
@inline function cons2entropy(u, aux_vars,
                              equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    return SVector(u[1], z, z)
end

# Flux as a function of the state vector u, as well as the auxiliary variables aux_vars, 
# which contain the geometric information required for the covariant form
@inline function flux(u, aux_vars, orientation::Integer,
                      equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    J = area_element(aux_vars, equations)
    vcon = velocity_contravariant(u, equations)
    return SVector(J * u[1] * vcon[orientation], z, z)
end

@inline function flux(u, aux_vars, normal_direction::AbstractVector,
                      equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    J = area_element(aux_vars, equations)
    vcon = velocity_contravariant(u, equations)
    a = dot(vcon, normal_direction) # velocity in normal direction
    return SVector(J * u[1] * a, z, z)
end

# Local Lax-Friedrichs dissipation which is not applied to the contravariant velocity 
# components, as they should remain unchanged in time
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr, aux_vars_ll,
                                                              aux_vars_rr,
                                                              orientation_or_normal_direction,
                                                              equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u_ll))
    J = area_element(aux_vars_ll, equations)
    λ = dissipation.max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                  orientation_or_normal_direction, equations)
    return -0.5f0 * J * λ * SVector(u_rr[1] - u_ll[1], z, z)
end

# Maximum contravariant wave speed with respect to specific basis vector
@inline function max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               orientation::Integer,
                               equations::CovariantLinearAdvectionEquation2D)
    vcon_ll = velocity_contravariant(u_ll, equations)  # Contravariant components on left side
    vcon_rr = velocity_contravariant(u_rr, equations)  # Contravariant components on right side
    return max(abs(vcon_ll[orientation]), abs(vcon_rr[orientation]))
end

@inline function max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                     normal_direction::AbstractVector,
                                     equations::CovariantLinearAdvectionEquation2D)
    vcon_ll = velocity_contravariant(u_ll, equations)  # Contravariant components on left side
    vcon_rr = velocity_contravariant(u_rr, equations)  # Contravariant components on right side
    # Calculate the velocity in the normal direction
    a_ll = abs(dot(vcon_ll, normal_direction))
    a_rr = abs(dot(vcon_rr, normal_direction))
    return max(a_ll, a_rr)
end

# Maximum wave speeds in each direction for CFL calculation
@inline function max_abs_speeds(u, aux_vars,
                                equations::CovariantLinearAdvectionEquation2D)
    return abs.(velocity_contravariant(u, equations))
end

# If the initial velocity field is defined in Cartesian coordinates and the chosen global 
# coordinate system is spherical, perform the appropriate conversion
@inline function cartesian2global(u, x,
                                  ::CovariantLinearAdvectionEquation2D{GlobalSphericalCoordinates})
    vlon, vlat, vrad = cartesian2spherical(u[2], u[3], u[4], x)
    return SVector(u[1], vlon, vlat, vrad)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is Cartesian, perform the appropriate conversion
@inline function spherical2global(u, x,
                                  ::CovariantLinearAdvectionEquation2D{GlobalCartesianCoordinates})
    vx, vy, vz = spherical2cartesian(u[2], u[3], u[4], x)
    return SVector(u[1], vx, vy, vz)
end

# If the initial velocity field is defined in spherical coordinates and the chosen global 
# coordinate system is spherical, do not convert
@inline function spherical2global(u, x,
                                  ::CovariantLinearAdvectionEquation2D{GlobalSphericalCoordinates})
    return u
end
end # @muladd
