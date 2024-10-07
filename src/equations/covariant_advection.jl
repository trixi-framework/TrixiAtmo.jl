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

# Directional flux that takes in the normal components in reference space
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            ::CovariantLinearAdvectionEquation2D,
                            elements, i, j, element)
    z = zero(eltype(u))
    v_n = u[2] * normal_direction[1] + u[3] * normal_direction[2]
    J = volume_element(elements, i, j, element)
    return SVector(J * u[1] * v_n, z, z)
end

# Local Lax-Friedrichs dissipation which is not applied to the contravariant velocity 
# components, as they should remain unchanged in time
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              normal_direction::AbstractVector,
                                                              equations::CovariantLinearAdvectionEquation2D,
                                                              elements, i, j, element)
    z = zero(eltype(u_ll))
    J = volume_element(elements, i, j, element)
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations,
                                  elements, i, j, element)
    return -0.5f0 * J * λ * SVector(u_rr[1] - u_ll[1], z, z)
end

# Maximum wave speed in the normal direction
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           ::CovariantLinearAdvectionEquation2D,
                                           elements, i, j, element)
    v_n_ll = u_ll[2] * normal_direction[1] + u_ll[3] * normal_direction[2]
    v_n_rr = u_rr[2] * normal_direction[1] + u_rr[3] * normal_direction[2]
    return max(abs(v_n_ll), abs(v_n_rr))
end

# Maximum wave speeds with respect to the contravariant basis
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
