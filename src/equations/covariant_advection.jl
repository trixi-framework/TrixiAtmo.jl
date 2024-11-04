@muladd begin
#! format: noindent

@doc raw"""
    CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations{2, 3, 3}

A variable-coefficient linear advection equation can be defined on a two-dimensional
manifold $S \subset \mathbb{R}^3$ as
```math
\partial_t h + \nabla_S \cdot (h \vec{v}) = 0.
```
We treat this problem as a system of equations in which the first variable is the scalar
conserved quantity $h$, and the second two are the contravariant components $v^1$ and $v^2$ 
used in the expansion 
```math
\vec{v} = v^1 \vec{a}_1 + v^2 \vec{a}_2,
```
which are spatially varying but remain constant in time (i.e. no flux or dissipation is 
applied to such variables). The resulting system is then given on the reference element as 
```math
J \frac{\partial}{\partial t}\left[\begin{array}{c} h \\ v^1 \\ v^2 \end{array}\right] +
\frac{\partial}{\partial \xi^1} \left[\begin{array}{c} J h v^1 \\ 0 \\ 0 \end{array}\right]
+ 
\frac{\partial}{\partial \xi^2} \left[\begin{array}{c} J h v^2 \\ 0 \\ 0 \end{array}\right] 
= \left[\begin{array}{c} 0 \\ 0 \\ 0 \end{array}\right],
```
where $J$ is the area element (see the documentation for [`P4estElementContainerCovariant`](@ref)).
!!! note
    The initial condition should be prescribed as $[h, u, v]^{\mathrm{T}}$ in terms of the 
    global velocity components $u$ and $v$ (i.e. zonal and meridional components in the 
    case of a spherical shell). The transformation to local contravariant components $v^1$ 
    and $v^2$ is handled internally within `Trixi.compute_coefficients!`.
"""
struct CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations{2, 3, 3} end

function Trixi.varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation2D)
    return ("scalar", "v_con_1", "v_con_2")
end

# Compute the entropy variables (requires element container and indices)
@inline function Trixi.cons2entropy(u, equations::CovariantLinearAdvectionEquation2D,
                                    elements, i, j, element)
    z = zero(eltype(u))
    return SVector{3}(u[1], z, z)
end

# The flux for the covariant form takes in the element container and node/element indices
# in order to give the flux access to the geometric information
@inline function Trixi.flux(u, orientation::Integer,
                            ::CovariantLinearAdvectionEquation2D,
                            elements, i, j, element)
    z = zero(eltype(u))
    J = volume_element(elements, i, j, element)
    return SVector(J * u[1] * u[orientation + 1], z, z)
end

# Local Lax-Friedrichs dissipation which is not applied to the contravariant velocity 
# components, as they should remain unchanged in time
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::CovariantLinearAdvectionEquation2D,
                                                              elements, i_ll, j_ll,
                                                              i_rr, j_rr, element)
    z = zero(eltype(u_ll))
    J = volume_element(elements, i_ll, j_ll, element)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations, elements, i_ll, j_ll, i_rr, j_rr, element)
    return -0.5f0 * J * λ * SVector(u_rr[1] - u_ll[1], z, z)
end

# Maximum wave speed with respect to the a specific orientation
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                           ::CovariantLinearAdvectionEquation2D,
                                           elements, i_ll, j_ll, i_rr, j_rr, element)
    return max(abs(u_ll[orientation + 1]), abs(u_rr[orientation + 1]))
end

# Maximum wave speeds in each direction for CFL calculation
@inline function Trixi.max_abs_speeds(u, ::CovariantLinearAdvectionEquation2D,
                                      elements, i, j, element)
    return abs(u[2]), abs(u[3])
end

# Convert contravariant velocity/momentum components to zonal and meridional components
@inline function contravariant2spherical(u::SVector{3},
                                         ::CovariantLinearAdvectionEquation2D,
                                         elements, i, j, element)
    v_lon, v_lat = contravariant2spherical(u[2], u[3], elements, i, j, element)
    return SVector(u[1], v_lon, v_lat)
end

# Convert zonal and meridional velocity/momentum components to contravariant components
@inline function spherical2contravariant(u::SVector{3},
                                         ::CovariantLinearAdvectionEquation2D,
                                         elements, i, j, element)
    v_con_1, v_con_2 = spherical2contravariant(u[2], u[3], elements, i, j, element)
    return SVector(u[1], v_con_1, v_con_2)
end
end # @muladd
