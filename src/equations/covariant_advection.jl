@muladd begin
#! format: noindent

"""
    Variable-coefficient linear advection equation in covariant form
"""
struct CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations2D{3} end

function Trixi.varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation2D)
    # The first variable is the scalar conserved quantity. 
    # The second two are the contravariant velocity components, 
    # which are spatially varying but remain constant in time.
    return ("scalar", "v_con_1", "v_con_2")
end

Trixi.cons2entropy(u, ::CovariantLinearAdvectionEquation2D) = u

@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              normal_direction::AbstractVector,
                                                              equations::CovariantLinearAdvectionEquation2D,
                                                              i, j, element, cache)
    z = zero(eltype(u_ll))
    J = 1 / cache.elements.inverse_jacobian[i, j, element]
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations, i, j,
                                  element, cache)
    return -0.5f0 * J * λ * SVector(u_rr[1] - u_ll[1], z, z)
end

@inline function Trixi.flux(u, orientation::Integer,
                            ::CovariantLinearAdvectionEquation2D,
                            i, j, element, cache)
    J = 1 / cache.elements.inverse_jacobian[i, j, element]
    z = zero(eltype(u))
    return SVector(J * u[orientation + 1] * u[1], z, z)
end

@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            ::CovariantLinearAdvectionEquation2D,
                            i, j, element, cache)
    z = zero(eltype(u))
    J = 1 / cache.elements.inverse_jacobian[i, j, element]
    v_n = u[2] * normal_direction[1] + u[3] * normal_direction[2]
    return SVector(J * v_n * u[1], z, z)
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           ::CovariantLinearAdvectionEquation2D,
                                           i, j, element, cache)
    v_n_ll = u_ll[2] * normal_direction[1] + u_ll[3] * normal_direction[2]
    v_n_rr = u_rr[2] * normal_direction[1] + u_rr[3] * normal_direction[2]
    return max(abs(v_n_ll), abs(v_n_rr))
end

# Maximum wave speeds with respect to the contravariant basis
@inline function Trixi.max_abs_speeds(u, ::CovariantLinearAdvectionEquation2D,
                                      i, j, element, cache)
    return abs(u[2]), abs(u[3])
end
end # @muladd
