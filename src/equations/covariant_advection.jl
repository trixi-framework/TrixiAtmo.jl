###############################################################################
# Linear advection equation in general covariant form
###############################################################################

struct CovariantLinearAdvectionEquation3D <: AbstractCovariantEquations3D{4} end

function Trixi.varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation3D)
    return ("contravariant vel. 1", "contravariant vel. 2", "contravariant vel. 3", 
        "density")
end

default_analysis_integrals(::CovariantLinearAdvectionEquation3D) = ()
default_analysis_errors(::CovariantLinearAdvectionEquation3D) = ()

"""
    Custom dissipation to ensure no flux is applied to variable coefficients
"""
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                        normal_direction::AbstractVector,
                        equations::CovariantLinearAdvectionEquation3D)
    z = zero(eltype(u_ll))
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations)
    return -0.5f0 * λ * SVector(z, z, z, u_rr[4] - u_ll[4])
end

"""
    Compute a given contravariant flux component
"""
@inline function Trixi.flux(u, orientation::Int, ::CovariantLinearAdvectionEquation3D)
    z = zero(eltype(u))
    return SVector(z, z, z, u[orientation] * u[4])
end

"""
    Compute the flux component in a given normal direction
"""
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                               ::CovariantLinearAdvectionEquation3D)
    z = zero(eltype(u))
    return SVector(z, z, z, 
                   (u[1] * normal_direction[1] + u[2] * normal_direction[2]) * u[4])
end

"""
    Maximum directional wave speed for Lax-Friedrichs dissipation
"""
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                           ::CovariantLinearAdvectionEquation3D)

    v_n_ll = u_ll[1] * normal_direction[1] + u_ll[2] * normal_direction[2]
    v_n_rr = u_rr[1] * normal_direction[1] + u_rr[2] * normal_direction[2]
    return max(abs(v_n_ll), abs(v_n_rr))
end

"""
    Maximum wave speeds with respect to the contravariant basis
"""
@inline function Trixi.max_abs_speeds(u, ::CovariantLinearAdvectionEquation3D)
    return abs(u[1]), abs(u[2]), zero(eltype(u))
end