###############################################################################
# Linear advection equation in covariant form
###############################################################################

@muladd begin
#! format: noindent

struct CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations2D{3} end

"""
    The first variable is the scalar conserved quantity. The second two are the contravariant velocity components, which are spatially varying but remain constant in time.
"""
function Trixi.varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation2D)
    return ("scalar", "v_con_1", "v_con_2")
end

Trixi.cons2entropy(u, ::CovariantLinearAdvectionEquation2D) = u

"""
    Custom dissipation to ensure no flux is applied to variable coefficients
"""
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                        normal_direction::AbstractVector,
                        equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u_ll))
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations)
    return -0.5f0 * λ * SVector(u_rr[1] - u_ll[1], z, z)
end

"""
    Compute a given contravariant flux component
"""
@inline function Trixi.flux(u, orientation::Int, ::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    return SVector(u[orientation+1] * u[1], z, z)
end

"""
    Compute the flux component in a given normal direction
"""
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                               ::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    return SVector((u[2] * normal_direction[1] + u[3] * normal_direction[2]) * u[1], z, z)
end

"""
    Maximum directional wave speed for Lax-Friedrichs dissipation
"""
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           ::CovariantLinearAdvectionEquation2D)
    v_n_ll = u_ll[2] * normal_direction[1] + u_ll[3] * normal_direction[2]
    v_n_rr = u_rr[2] * normal_direction[1] + u_rr[3] * normal_direction[2]
    return max(abs(v_n_ll), abs(v_n_rr))
end

"""
    Maximum wave speeds with respect to the contravariant basis
"""
@inline function Trixi.max_abs_speeds(u, ::CovariantLinearAdvectionEquation2D)
    return abs(u[2]), abs(u[3])
end
end # @muladd