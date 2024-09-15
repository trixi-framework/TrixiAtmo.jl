@muladd begin
#! format: noindent

# Variable-coefficient linear advection equation in covariant form on a two-dimensional
# surface in three-dimensional space
struct CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations{2, 3, 3} end

# The first variable is the scalar conserved quantity, and the second two are the 
# contravariant velocity components, which are spatially varying but constant in time.
function Trixi.varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation2D)
    return ("scalar", "v_con_1", "v_con_2")
end

# We will measure the entropy/energy simply as half of the square of the scalar variable
Trixi.cons2entropy(u, ::CovariantLinearAdvectionEquation2D) = SVector{3}(u[1],
                                                                         zero(eltype(u)),
                                                                         zero(eltype(u)))

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
