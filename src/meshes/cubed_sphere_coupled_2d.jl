###############################################################################
# Geometry and connectivity information for the cubed sphere
###############################################################################

"""
    A single face of the cubed sphere (Sadourny, 1972) of radius 
    a ∈ R⁺ embedded in R³, where the mapping from parameter space 
    x¹,x² ∈ [-π/4,π/4] employs an equiangular gnomonic projection 
    (Ronchi et al., 1996)
"""
struct CubedSphereFace2D{face_number, RealT <: Real} <: AbstractManifold{2,3} 
    a::RealT  # radius of the sphere
    function CubedSphereFace2D{face_number}(
        a::RealT) where {face_number, RealT <: Real}
        return new{face_number,RealT}(a)
    end
end

"""
    Connectivity for the cubed sphere as a vector of six NamedTuples of 
    BoundaryConditionCoupled structs
"""
function cubed_sphere_connectivity(coupling_function)
    return [(  # face 1 (+x)
        x_neg = BoundaryConditionCoupled(
            4, (:end, :i_forward), Float64, coupling_function),
        x_pos = BoundaryConditionCoupled(
            2, (:begin, :i_forward), Float64, coupling_function),
        y_neg = BoundaryConditionCoupled(
            6, (:i_forward, :end), Float64, coupling_function),
        y_pos = BoundaryConditionCoupled(
            5, (:i_forward, :begin), Float64, coupling_function)),
    (   # face 2 (+y)
        x_neg = BoundaryConditionCoupled(
            1, (:end, :i_forward), Float64, coupling_function),
        x_pos = BoundaryConditionCoupled(
            3, (:begin, :i_forward), Float64, coupling_function),
        y_neg = BoundaryConditionCoupled(
            6, (:end, :i_backward), Float64, coupling_function),
        y_pos = BoundaryConditionCoupled(
            5, (:end,:i_forward), Float64, coupling_function)),
    (   # face 3 (-x)
        x_neg = BoundaryConditionCoupled(
            2, (:end, :i_forward), Float64, coupling_function),
        x_pos = BoundaryConditionCoupled(
            4, (:begin, :i_forward), Float64, coupling_function),
        y_neg = BoundaryConditionCoupled(
            6, (:i_backward, :begin), Float64, coupling_function),
        y_pos = BoundaryConditionCoupled(
            5, (:i_backward, :end), Float64, coupling_function)),
    (   # face 4 (-y)
        x_neg = BoundaryConditionCoupled(
            3, (:end, :i_forward), Float64, coupling_function),
        x_pos = BoundaryConditionCoupled(
            1, (:begin, :i_forward), Float64, coupling_function),
        y_neg = BoundaryConditionCoupled(
            6, (:begin, :i_forward), Float64, coupling_function),
        y_pos= BoundaryConditionCoupled(
            5, (:begin,:i_backward), Float64, coupling_function)),
    (  # face 5 (+z)
        x_neg = BoundaryConditionCoupled(
            4, (:i_backward, :end), Float64, coupling_function),
        x_pos = BoundaryConditionCoupled(
            2, (:i_forward, :end), Float64, coupling_function),
        y_neg = BoundaryConditionCoupled(
            1, (:i_forward, :end), Float64, coupling_function),
        y_pos = BoundaryConditionCoupled(
            3, (:i_backward, :end), Float64, coupling_function)),
    (   # face 6 (-z)
        x_neg = BoundaryConditionCoupled(
            4, (:i_forward, :begin), Float64, coupling_function),
        x_pos = BoundaryConditionCoupled(
            2, (:i_backward, :begin), Float64, coupling_function),
        y_neg = BoundaryConditionCoupled(
            3, (:i_backward, :begin), Float64, coupling_function),
        y_pos = BoundaryConditionCoupled(
            1, (:i_forward, :begin), Float64, coupling_function))]
end

@inline basis_covariant(x, manifold) = inv(basis_contravariant(x, manifold))
@inline function metric_covariant(x, manifold::CubedSphereFace2D)
    A = basis_covariant(x, manifold)
    return A' * A
end

@inline metric_contravariant(x, manifold) = inv(metric_covariant(x, manifold))

@inline function area_element(x, manifold::CubedSphereFace2D)
    δ = sqrt(1 + tan(x[1])^2 + tan(x[2])^2)
    return manifold.a^2 / (δ^3 *(cos(x[1])*cos(x[2]))^2)
end

"""
     Transform longitude λ and latitude θ to global Cartesian
    coordinates x, y,z
"""
@inline function global2cartesian(λ, θ, manifold::CubedSphereFace2D)
    (; a) = manifold
    return a*cos(λ)*cos(θ), a*sin(λ)*cos(θ), a*sin(θ)
end

"""
     Transform global Cartesian coordinates x, y, z to longitude λ
    and latitude θ
"""
@inline function cartesian2global(x, y, z, manifold::CubedSphereFace2D)
    return atan(y,x), asin(z/manifold.a)
end

"""
     Compute the Christoffel Symbols of the second kind
"""
@inline function christoffel_symbols(x, ::CubedSphereFace2D)

    X, Y = tan(x[1]), tan(x[2])
    scl = 1/(1 + X^2 + Y^2)

    Γ1_11 = 2*scl*X*Y^2
    Γ1_12 = -scl*Y*(1+Y^2)
    Γ1_21 = Γ1_12
    Γ1_22 = 0.0f0

    Γ2_11 = 0.0f0
    Γ2_12 = -scl*X*(1+X^2)
    Γ2_21 = Γ2_12
    Γ2_22 = 2*scl*X^2*Y
    
    return SMatrix{2,2}(Γ1_11, Γ1_21, Γ1_12, Γ1_22), 
           SMatrix{2,2}(Γ2_11, Γ2_21, Γ2_12, Γ2_22)
end

###############################################################################
# coordinate transformations for face 1 (+x) 
###############################################################################
@inline function face2cartesian(x, manifold::CubedSphereFace2D{1})
    return global2cartesian(x[1], atan(tan(x[2]),sec(x[1])), manifold)
end
@inline function face2global(x, ::CubedSphereFace2D{1}) 
    return x[1], atan(tan(x[2]),sec(x[1]))
end
@inline function global2face(λ, θ, ::CubedSphereFace2D{1})
    return atan(tan(λ)), atan(tan(θ)*sec(λ))
end
@inline function basis_contravariant(x, manifold::CubedSphereFace2D{1})
    λ_loc, θ_loc = x[1], atan(tan(x[2]), sec(x[1]))
    cos²x1, cos²x2 = cos(x[1])^2, cos(x[2])^2
    secλ_loc, secθ_loc = sec(λ_loc), sec(θ_loc)
    return secθ_loc*secλ_loc/manifold.a * SMatrix{2,2}( 
        cos²x1 * secλ_loc, cos²x2*tan(θ_loc)*tan(λ_loc), 
        0, cos²x2*secθ_loc)
end

###############################################################################
# coordinate transformations for face 2 (+y)
###############################################################################
@inline function face2cartesian(x, manifold::CubedSphereFace2D{2})
    return global2cartesian(x[1] + π/2, atan(tan(x[2]),sec(x[1])), manifold)
end
@inline function face2global(x, ::CubedSphereFace2D{2}) 
    return x[1] + π/2, atan(tan(x[2]),sec(x[1]))
end
@inline function global2face(λ, θ, ::CubedSphereFace2D{2})
    λ_loc = λ - π/2
    return atan(tan(λ_loc)), atan(tan(θ)*sec(λ_loc))
end
@inline function basis_contravariant(x, manifold::CubedSphereFace2D{2})    
    λ_loc, θ_loc = x[1], atan(tan(x[2]), sec(x[1]))
    cos²x1, cos²x2 = cos(x[1])^2, cos(x[2])^2
    secλ_loc, secθ_loc = sec(λ_loc), sec(θ_loc)
    return secθ_loc*secλ_loc/manifold.a * SMatrix{2,2}( 
        cos²x1 * secλ_loc, cos²x2*tan(θ_loc)*tan(λ_loc), 
        0, cos²x2*secθ_loc)
end

###############################################################################
# coordinate transformations for face 3 (-x)
###############################################################################
@inline function face2cartesian(x, manifold::CubedSphereFace2D{3})
    return global2cartesian(x[1] + π, atan(tan(x[2]),sec(x[1])), manifold)
end
@inline function face2global(x, ::CubedSphereFace2D{3}) 
    return x[1] + π, atan(tan(x[2]),sec(x[1]))
end
@inline function global2face(λ, θ, ::CubedSphereFace2D{3})
    λ_loc = λ - π
    return atan(tan(λ_loc)), atan(tan(θ)*sec(λ_loc))
end
@inline function basis_contravariant(x, manifold::CubedSphereFace2D{3})  
    λ_loc, θ_loc = x[1], atan(tan(x[2]), sec(x[1]))
    cos²x1, cos²x2 = cos(x[1])^2, cos(x[2])^2
    secλ_loc, secθ_loc = sec(λ_loc), sec(θ_loc)
    return secθ_loc*secλ_loc/manifold.a * SMatrix{2,2}( 
        cos²x1 * secλ_loc, cos²x2*tan(θ_loc)*tan(λ_loc), 
        0, cos²x2*secθ_loc)
end

###############################################################################
# coordinate transformations for face 4 (-y)
###############################################################################
@inline function face2cartesian(x, manifold::CubedSphereFace2D{4})
    return global2cartesian(x[1] + 3π/2, 
        atan(tan(x[2]),sec(x[1])), manifold)
end
@inline function face2global(x, ::CubedSphereFace2D{4}) 
    return x[1] + 3π/2, atan(tan(x[2]),sec(x[1]))
end
@inline function global2face(λ, θ, ::CubedSphereFace2D{4})
    λ_loc = λ - 3π/2
    return atan(tan(λ_loc)), atan(tan(θ)*sec(λ_loc))
end
@inline function basis_contravariant(x, manifold::CubedSphereFace2D{4})    
    λ_loc, θ_loc = x[1], atan(tan(x[2]), sec(x[1]))
    cos²x1, cos²x2 = cos(x[1])^2, cos(x[2])^2
    secλ_loc, secθ_loc = sec(λ_loc), sec(θ_loc)
    return secθ_loc*secλ_loc/manifold.a * SMatrix{2,2}( 
        cos²x1 * secλ_loc, cos²x2*tan(θ_loc)*tan(λ_loc), 
        0, cos²x2*secθ_loc)
end

###############################################################################
# coordinate transformations for face 5 (+z)
###############################################################################
@inline function face2cartesian(x, manifold::CubedSphereFace2D{5})
    x_c, y_c, z_c = global2cartesian(x[1], 
        atan(tan(x[2]),sec(x[1])), manifold)
    return -z_c, y_c, x_c
end
@inline function face2global(x, manifold::CubedSphereFace2D{5}) 
    θ_loc = atan(tan(x[2]),sec(x[1]))
    x_loc, y_loc, z_loc = global2cartesian(x[1], θ_loc, manifold)
    return cartesian2global(-z_loc, y_loc, x_loc, manifold)
end
@inline function global2face(λ, θ, ::CubedSphereFace2D{5})
    return atan(sin(λ)*cot(θ)), atan(-cos(λ)*cot(θ))
end
@inline function basis_contravariant(x, manifold::CubedSphereFace2D{5})  
    λ, θ = face2global(x, manifold)
    cos²x1, cos²x2 = cos(x[1])^2, cos(x[2])^2
    sinλ, sinθ, cosλ = sin(λ), sin(θ), cos(λ)
    return 1/(manifold.a*sinθ^2) * SMatrix{2,2}(
        cos²x1*sinθ*cosλ, cos²x2*sinθ*sinλ, -cos²x1*sinλ, cos²x2*cosλ)
end

###############################################################################
# coordinate transformations for face 6 (-z)
###############################################################################
@inline function face2cartesian(x, manifold::CubedSphereFace2D{6})
    x_c, y_c, z_c = global2cartesian(x[1], 
        atan(tan(x[2]),sec(x[1])), manifold)
    return z_c, y_c, -x_c
end
@inline function face2global(x, manifold::CubedSphereFace2D{6}) 
    θ_loc = atan(tan(x[2]),sec(x[1]))
    x_loc, y_loc, z_loc = global2cartesian(x[1], θ_loc, manifold)
    return cartesian2global(z_loc, y_loc, -x_loc, manifold)
end
@inline function global2face(λ, θ, ::CubedSphereFace2D{6})
    return atan(-sin(λ)*cot(θ)), atan(-cos(λ)*cot(θ)) 
end
@inline function basis_contravariant(x, manifold::CubedSphereFace2D{6})    
    λ, θ = face2global(x, manifold)
    cos²x1, cos²x2 = cos(x[1])^2, cos(x[2])^2
    sinλ, sinθ, cosλ = sin(λ), sin(θ), cos(λ)
    return 1/(manifold.a*sinθ^2) * SMatrix{2,2}(
        -cos²x1*sinθ*cosλ, cos²x2*sinθ*sinλ, cos²x1*sinλ, cos²x2*cosλ)
end