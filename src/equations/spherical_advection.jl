###############################################################################
# Linear advection equation on the sphere
###############################################################################

struct SphericalLinearAdvectionEquation2D{Manifold <: 
    AbstractManifold{2,3}} <: Trixi.AbstractEquations{2,4}
    manifold::Manifold
end

Trixi.varnames(::typeof(cons2cons), ::SphericalLinearAdvectionEquation2D) = (
    "H",  # height scaled by area element
    "sqrtG",  # area element
    "v1", "v2")  # contravariant velocity components

Trixi.varnames(::typeof(cons2prim), ::SphericalLinearAdvectionEquation2D) = (
    "h",  # height
    "sqrtG",  # area 
    "v1", "v2")  # contravariant velocity components

default_analysis_integrals(::SphericalLinearAdvectionEquation2D) = ()
default_analysis_errors(::SphericalLinearAdvectionEquation2D) = ()

"""
    Convenience functions to extract named variables
"""
@inline scaled_waterheight(u, ::SphericalLinearAdvectionEquation2D) = u[1]
@inline area_element(u, ::SphericalLinearAdvectionEquation2D) = u[2]
@inline waterheight(u, ::SphericalLinearAdvectionEquation2D) = u[1]/u[2]
@inline velocity_contravariant(u, ::SphericalLinearAdvectionEquation2D) = 
    SVector(u[3], u[4])

""" 
    Convert conservative to primitive variables
"""
@inline Trixi.cons2prim(u, ::SphericalLinearAdvectionEquation2D) = 
    SVector(u[1]/u[2], u[2], u[3], u[4])

""" 
    Convert primitive to conservative variables
"""
@inline Trixi.prim2cons(q, ::SphericalLinearAdvectionEquation2D) = 
    SVector(u[1]*u[2], u[2], u[3], u[4])

"""
    Convert global variables (i.e. velocities in longitude/latitude
    coordinates) to conservative variables
"""
@inline function global2cons(h, v_λ, v_θ, λ, θ,
    equations::SphericalLinearAdvectionEquation2D)

    # Get local coordinates on a given face
    x = global2face(λ, θ, equations.manifold)

    # Get covariant basis matrix
    A = basis_covariant(x, equations.manifold)

    # Compute contravariant velocity components
    v1, v2 = A \ SVector(v_λ, v_θ)

    # Compute the area element
    G = A' * A
    sqrtG = sqrt(G[1,1]*G[2,2] - G[1,2]^2)

    # return the conservative variables
    return SVector(sqrtG * h, sqrtG, v1, v2)
end

"""
    Coupling function for multiple manifolds. Transform the velocity vector from
    the "other" coordinate system to the "own" coordinate system 
"""
@inline function spherical_coupling(x_other, u_other,
    equations_other::SphericalLinearAdvectionEquation2D, 
    equations_own::SphericalLinearAdvectionEquation2D)
    
    # get the water height on "other" face
    h = waterheight(u_other, equations_other)

    # Get global latitude-longitude coordinates
    λ, θ = face2global(x_other, equations_other.manifold)

    # Get covariant basis matrix on "other" face
    A = basis_covariant(x_other, equations_other.manifold)

    # transform to global zonal and meridional velocity components
    v_λ, v_θ = A * velocity_contravariant(u_other, equations_other)

    # convert from global coordinates to local coordinates on own manifold
    return global2cons(h, v_λ, v_θ, λ, θ, equations_own)
end

"""
    Custom dissipation to ensure no flux is applied to variable coefficients
"""
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
    normal_direction::AbstractVector,
    equations::SphericalLinearAdvectionEquation2D)

    z = zero(eltype(u_ll))
    λ = dissipation.max_abs_speed(u_ll, u_rr, normal_direction, equations)
    return -0.5f0 * λ * SVector(u_rr[1] - u_ll[1], z, z, z)
end

"""
    Compute a given contravariant flux component
"""
@inline function Trixi.flux(u, orientation::Int, 
    equations::SphericalLinearAdvectionEquation2D)

    z = zero(eltype(u))
    v_con = velocity_contravariant(u, equations)

    return SVector(v_con[orientation] * scaled_waterheight(u,equations), z, z, z)
end

"""
    Compute the flux component in a given normal direction
"""
@inline function Trixi.flux(u, normal_direction::AbstractVector,
    equations::SphericalLinearAdvectionEquation2D)

    z = zero(eltype(u))
    v_n = dot(velocity_contravariant(u, equations), normal_direction)

    return SVector(v_n * scaled_waterheight(u,equations), z, z, z)
end

"""
    Maximum directional wave speed for Lax-Friedrichs dissipation
"""
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction,
    equations::SphericalLinearAdvectionEquation2D)

    v_n_ll = dot(velocity_contravariant(u_ll, equations), normal_direction)
    v_n_rr = dot(velocity_contravariant(u_rr, equations), normal_direction)
    
    return max(abs(v_n_ll), abs(v_n_rr))
end

"""
    Maximum wave speeds with respect to the contravariant basis 
    (used in CFL-based time step control)
"""
@inline function Trixi.max_abs_speeds(u,
    equations::SphericalLinearAdvectionEquation2D)
    
    v = velocity_contravariant(u, equations)

    return abs(v[1]), abs(v[2])
end
