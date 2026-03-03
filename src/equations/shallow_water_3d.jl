@muladd begin
#! format: noindent

@doc raw"""
    ShallowWaterEquations3D(; gravity, rotation_rate = 0, H0 = 0)

Rotating shallow water equations (SWE) in three space dimensions to be solved on curved manifolds (e.g., on a spherical shell). 
The equations are given by
```math
\begin{aligned}
  \frac{\partial h}{\partial t} + \frac{\partial}{\partial x}(h v_1)
    + \frac{\partial}{\partial y}(h v_2) + \frac{\partial}{\partial z}(h v_3) &= s_h \\
    \frac{\partial}{\partial t}(h v_1) + \frac{\partial}{\partial x}\left(h v_1^2 + \frac{g}{2}h^2\right)
    + \frac{\partial}{\partial y}(h v_1 v_2) + \frac{\partial}{\partial z}(h v_1 v_3) g h \frac{\partial b}{\partial x} &= s_{hv_1} \\
    \frac{\partial}{\partial t}(h v_2) + \frac{\partial}{\partial x}(h v_1 v_2)
    + \frac{\partial}{\partial y}\left(h v_2^2 + \frac{g}{2}h^2\right) + \frac{\partial}{\partial z}(h v_2 v_3) + g h \frac{\partial b}{\partial y} &= s_{hv_2} \\
    \frac{\partial}{\partial t}(h v_3) + \frac{\partial}{\partial x}(h v_1 v_3)
    + \frac{\partial}{\partial y}(h v_2 v_3) + \frac{\partial}{\partial z}\left(h v_3^2 + \frac{g}{2}h^2\right) + g h \frac{\partial b}{\partial z}&= s_{hv_3}.
\end{aligned}
```
The unknown quantities of the SWE are the water height ``h`` and the velocities ``\mathbf{v} = (v_1, v_2, v_3)^T``.
The gravitational acceleration is denoted by `g`.

The 3D Shallow Water Equations (SWE) extend the 2D SWE to model shallow water flows on 2D manifolds embedded within 3D space. 
To confine the flow to the 2D manifold, a source term incorporating a Lagrange multiplier is applied to the momentum equations using the function [`source_terms_lagrange_multiplier`](@ref). 
This term effectively removes momentum components that are normal to the manifold, ensuring the flow remains 
constrained within the 2D surface.

To incorporate the effect of the rotation of the manifold, use the function [`source_terms_coriolis`](@ref), which adds the necessary Coriolis source terms to the momentum equations assuming a rotation around the ``z`` axis with a rotation rate in radians per time unit given by `rotation_rate`. To incorporate both Coriolis forces and the Lagrange multiplier terms, use [`source_terms_coriolis_lagrange_multiplier`](@ref).

The additional quantity ``H_0`` is also available to store a reference value for the total water height that
is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

In addition to the unknowns, TrixiAtmo.jl currently stores the bottom topography values at the approximation points
despite being fixed in time. This is done for convenience of computing the bottom topography gradients
on the fly during the approximation as well as computing auxiliary quantities like the total water height ``H``
or the entropy variables.
This affects the implementation and use of these equations in various ways:
* The flux values corresponding to the bottom topography must be zero.
* The bottom topography values must be included when defining initial conditions, boundary conditions or
  source terms.
* [`AnalysisCallback`](https://trixi-framework.github.io/TrixiDocumentation/stable/reference-trixi/#Trixi.AnalysisCallback) analyzes this variable.
* Trixi.jl's visualization tools will visualize the bottom topography by default.

References:
- J. Coté (1988). "A Lagrange multiplier approach for the metric terms of semi-Lagrangian models on the sphere". 
  Quarterly Journal of the Royal Meteorological Society 114, 1347-1352. [DOI: 10.1002/qj.49711448310](https://doi.org/10.1002/qj.49711448310)
- F. X. Giraldo (2001). "A spectral element shallow water model on spherical geodesic grids". 
  [DOI: 10.1002/1097-0363(20010430)35:8<869::AID-FLD116>3.0.CO;2-S](https://doi.org/10.1002/1097-0363(20010430)35:8%3C869::AID-FLD116%3E3.0.CO;2-S)
"""
struct ShallowWaterEquations3D{RealT <: Real} <:
       AbstractShallowWaterEquations{3, 5}
    gravity::RealT # gravitational acceleration
    rotation_rate::RealT  # rotation rate around z axis for Coriolis term 
    H0::RealT      # constant "lake-at-rest" total water height
end

# Allow for flexibility to set the gravitational acceleration within an elixir depending on the
# application where `gravity=1.0` or `gravity=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases.
function ShallowWaterEquations3D(; gravity, rotation_rate = zero(gravity),
                                 H0 = zero(gravity))
    ShallowWaterEquations3D(gravity, rotation_rate, H0)
end

have_nonconservative_terms(::ShallowWaterEquations3D) = True()
varnames(::typeof(cons2cons), ::ShallowWaterEquations3D) = ("h", "h_v1", "h_v2", "h_v3",
                                                            "b")
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
varnames(::typeof(cons2prim), ::ShallowWaterEquations3D) = ("H", "v1", "v2", "v3", "b")
# Calculate 1D flux for a single point
# Note, the bottom topography has no flux
@inline function flux(u, orientation::Integer, equations::ShallowWaterEquations3D)
    h, h_v1, h_v2, h_v3, _ = u
    v1, v2, v3 = velocity(u, equations)

    p = 0.5f0 * equations.gravity * h^2
    if orientation == 1
        f1 = h_v1
        f2 = h_v1 * v1 + p
        f3 = h_v1 * v2
        f4 = h_v1 * v3
    elseif orientation == 2
        f1 = h_v2
        f2 = h_v2 * v1
        f3 = h_v2 * v2 + p
        f4 = h_v2 * v3
    else # orientation == 3
        f1 = h_v3
        f2 = h_v3 * v1
        f3 = h_v3 * v2
        f4 = h_v3 * v3 + p
    end
    return SVector(f1, f2, f3, f4, zero(eltype(u)))
end

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized and the bottom topography has no flux
@inline function flux(u, normal_direction::AbstractVector,
                      equations::ShallowWaterEquations3D)
    h = waterheight(u, equations)
    v1, v2, v3 = velocity(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
               v3 * normal_direction[3]
    h_v_normal = h * v_normal
    p = 0.5f0 * equations.gravity * h^2

    f1 = h_v_normal
    f2 = h_v_normal * v1 + p * normal_direction[1]
    f3 = h_v_normal * v2 + p * normal_direction[2]
    f4 = h_v_normal * v3 + p * normal_direction[3]
    return SVector(f1, f2, f3, f4, zero(eltype(u)))
end

"""
    flux_wintermeyer_etal(u_ll, u_rr,
                          normal_direction::AbstractVector,
                          equations::ShallowWaterEquations3D)

Total energy conservative (mathematical entropy for shallow water equations) split form.
When the bottom topography is nonzero this scheme will be well-balanced when used as a `volume_flux`.
For the `surface_flux` either [`flux_wintermeyer_etal`](@ref) or [`flux_fjordholm_etal`](@ref) can
be used to ensure well-balancedness and entropy conservation.

Further details are available in Theorem 1 of the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
"""
@inline function flux_wintermeyer_etal(u_ll, u_rr,
                                       normal_direction::AbstractVector,
                                       equations::ShallowWaterEquations3D)
    # Unpack left and right state
    h_ll, h_v1_ll, h_v2_ll, h_v3_ll, _ = u_ll
    h_rr, h_v1_rr, h_v2_rr, h_v3_rr, _ = u_rr

    # Get the velocities on either side
    v1_ll, v2_ll, v3_ll = velocity(u_ll, equations)
    v1_rr, v2_rr, v3_rr = velocity(u_rr, equations)

    # Average each factor of products in flux
    h_v1_avg = 0.5f0 * (h_v1_ll + h_v1_rr)
    h_v2_avg = 0.5f0 * (h_v2_ll + h_v2_rr)
    h_v3_avg = 0.5f0 * (h_v3_ll + h_v3_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * equations.gravity * h_ll * h_rr

    # Calculate fluxes depending on normal_direction
    f1 = h_v1_avg * normal_direction[1] + h_v2_avg * normal_direction[2] +
         h_v3_avg * normal_direction[3]
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

"""
    flux_fjordholm_etal(u_ll, u_rr,
                        normal_direction::AbstractVector,
                        equations::ShallowWaterEquations3D)

Total energy conservative (mathematical entropy for shallow water equations). When the bottom topography
is nonzero this should only be used as a surface flux otherwise the scheme will not be well-balanced.
For well-balancedness in the volume flux use [`flux_wintermeyer_etal`](@ref).

Details are available in Eq. (4.1) in the paper:
- Ulrik S. Fjordholm, Siddhartha Mishra and Eitan Tadmor (2011)
  Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography
  [DOI: 10.1016/j.jcp.2011.03.042](https://doi.org/10.1016/j.jcp.2011.03.042)
"""
@inline function flux_fjordholm_etal(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::ShallowWaterEquations3D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll, v3_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr, v3_rr = velocity(u_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]

    # Average each factor of products in flux
    h_avg = 0.5f0 * (h_ll + h_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    h2_avg = 0.5f0 * (h_ll^2 + h_rr^2)
    p_avg = 0.5f0 * equations.gravity * h2_avg
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)

    # Calculate fluxes depending on normal_direction
    f1 = h_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

"""
    flux_nonconservative_wintermeyer_etal(u_ll, u_rr, orientation::Integer,
                                          equations::ShallowWaterEquations3D)
    flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                          normal_direction::AbstractVector,
                                          equations::ShallowWaterEquations3D)

Non-symmetric two-point volume flux discretizing the nonconservative (source) term
that contains the gradient of the bottom topography [`ShallowWaterEquations3D`](@ref).

For the `surface_flux` either [`flux_wintermeyer_etal`](@ref) or [`flux_fjordholm_etal`](@ref) can
be used to ensure well-balancedness and entropy conservation.

Further details are available in the papers:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on
  curvilinear meshes
  [DOI: 10.48550/arXiv.2306.12699](https://doi.org/10.48550/arXiv.2306.12699)
"""
@inline function flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::ShallowWaterEquations3D)
    # Pull the necessary left and right state information
    h_ll = waterheight(u_ll, equations)
    b_jump = u_rr[5] - u_ll[5]

    # Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
    return SVector(0,
                   normal_direction[1] * equations.gravity * h_ll * b_jump,
                   normal_direction[2] * equations.gravity * h_ll * b_jump,
                   normal_direction[3] * equations.gravity * h_ll * b_jump,
                   0)
end

"""
    flux_nonconservative_fjordholm_etal(u_ll, u_rr, orientation::Integer,
                                        equations::ShallowWaterEquations3D)
    flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                        normal_direction::AbstractVector,
                                        equations::ShallowWaterEquations3D)

Non-symmetric two-point surface flux discretizing the nonconservative (source) term of
that contains the gradient of the bottom topography [`ShallowWaterEquations3D`](@ref).

This flux can be used together with [`flux_fjordholm_etal`](@ref) at interfaces to ensure entropy
conservation and well-balancedness.

Further details for the original finite volume formulation are available in
- Ulrik S. Fjordholm, Siddhartha Mishra and Eitan Tadmor (2011)
  Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography
  [DOI: 10.1016/j.jcp.2011.03.042](https://doi.org/10.1016/j.jcp.2011.03.042)
and for curvilinear 2D case in the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
"""
@inline function flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                     normal_direction::AbstractVector,
                                                     equations::ShallowWaterEquations3D)
    # Pull the necessary left and right state information
    h_ll, _, _, _, b_ll = u_ll
    h_rr, _, _, _, b_rr = u_rr

    h_average = 0.5f0 * (h_ll + h_rr)
    b_jump = b_rr - b_ll

    # Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
    f2 = normal_direction[1] * equations.gravity * h_average * b_jump
    f3 = normal_direction[2] * equations.gravity * h_average * b_jump
    f4 = normal_direction[3] * equations.gravity * h_average * b_jump

    # First and last equations do not have a nonconservative flux
    f1 = f5 = 0

    return SVector(f1, f2, f3, f4, f5)
end

"""
    source_terms_coriolis(u, du, x, t, 
                          equations::ShallowWaterEquations3D,
                          normal_direction)

Source term function to apply the Coriolis force with an angular velocity of
`equations.rotation_rate` around the ``z`` axis.

The vector `normal_direction` is perpendicular to the 2D manifold. By default, 
it is the normal contravariant basis vector.
"""
function source_terms_coriolis(u, du, x, t,
                               equations::ShallowWaterEquations3D,
                               normal_direction)
    _, h_v1, h_v2, h_v3, _ = u

    r2 = sum(normal_direction .^ 2) # Square of normal_direction's norm

    # Coriolis parameter scaled by norm of the normal vector
    # (projection of the angular velocity onto the normal vector to the manifold)
    f = 2 * equations.rotation_rate * normal_direction[3] / r2  # 2Ωsinθ / norm(normal_direction)

    # Compute the Coriolis source terms for the momentum equation as
    # s_mom = -f * cross(normal_direction, momentum)
    s2 = -f * (normal_direction[2] * h_v3 - normal_direction[3] * h_v2)
    s3 = -f * (normal_direction[3] * h_v1 - normal_direction[1] * h_v3)
    s4 = -f * (normal_direction[1] * h_v2 - normal_direction[2] * h_v1)

    return SVector(zero(eltype(u)), s2, s3, s4, zero(eltype(u)))
end

"""
    source_terms_lagrange_multiplier(u, du, x, t,
                                     equations::ShallowWaterEquations3D,
                                     normal_direction)

Source term function to apply a Lagrange multiplier to the semi-discretization
in order to constrain the momentum to a 2D manifold.

The vector `normal_direction` is perpendicular to the 2D manifold. By default, 
it is the normal contravariant basis vector.
"""
function source_terms_lagrange_multiplier(u, du, x, t,
                                          equations::ShallowWaterEquations3D,
                                          normal_direction)
    x_dot_div_f = (normal_direction[1] * du[2] +
                   normal_direction[2] * du[3] +
                   normal_direction[3] * du[4]) /
                  sum(normal_direction .^ 2)

    s2 = -normal_direction[1] * x_dot_div_f
    s3 = -normal_direction[2] * x_dot_div_f
    s4 = -normal_direction[3] * x_dot_div_f

    return SVector(0, s2, s3, s4, 0)
end

"""
    source_terms_coriolis_lagrange_multiplier(u, du, x, t,
                                              equations::ShallowWaterEquations3D,
                                              normal_direction)

Computes the Coriolis source term ([`source_terms_coriolis`](@ref)) and the
Lagrange multiplier source term ([`source_terms_lagrange_multiplier`](@ref)).
"""
function source_terms_coriolis_lagrange_multiplier(u, du, x, t,
                                                   equations::ShallowWaterEquations3D,
                                                   normal_direction)
    s_coriolis = source_terms_coriolis(u, du, x, t, equations, normal_direction)
    s_lagrange = source_terms_lagrange_multiplier(u, du, x, t, equations,
                                                  normal_direction)

    return s_coriolis + s_lagrange
end

"""
         clean_solution_lagrange_multiplier!(u, equations::ShallowWaterEquations3D, normal_direction)

Function to apply Lagrange multiplier discretely to the solution in order to constrain 
the momentum to a 2D manifold.

The vector `normal_direction` is perpendicular to the 2D manifold. By default, 
it is the normal contravariant basis vector.
"""
function clean_solution_lagrange_multiplier!(u, equations::ShallowWaterEquations3D,
                                             normal_direction)
    x_dot_div_f = (normal_direction[1] * u[2] +
                   normal_direction[2] * u[3] +
                   normal_direction[3] * u[4]) /
                  sum(normal_direction .^ 2)

    u[2] -= normal_direction[1] * x_dot_div_f
    u[3] -= normal_direction[2] * x_dot_div_f
    u[4] -= normal_direction[3] * x_dot_div_f
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::ShallowWaterEquations3D)
    # Extract and compute the velocities in the normal direction
    v1_ll, v2_ll, v3_ll = velocity(u_ll, equations)
    v1_rr, v2_rr, v3_rr = velocity(u_rr, equations)
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    # Compute the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)

    c_ll = sqrt(max(equations.gravity * h_ll, 0.0f0))
    c_rr = sqrt(max(equations.gravity * h_rr, 0.0f0))

    # The normal velocities are already scaled by the norm
    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::ShallowWaterEquations3D)
    # Extract and compute the velocities in the normal direction
    v1_ll, v2_ll, v3_ll = velocity(u_ll, equations)
    v1_rr, v2_rr, v3_rr = velocity(u_rr, equations)
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    # Compute the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)

    c_ll = sqrt(max(equations.gravity * h_ll, 0.0f0))
    c_rr = sqrt(max(equations.gravity * h_rr, 0.0f0))

    # The normal velocities are already scaled by the norm
    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_, abs(v_rr) + c_rr * norm_)
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the bottom topography
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::ShallowWaterEquations3D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations)
    diss = -0.5f0 * λ * (u_rr - u_ll)
    return SVector(diss[1], diss[2], diss[3], diss[4], zero(eltype(u_ll)))
end

# Specialization of [`DissipationLaxFriedrichsEntropyVariables`](https://trixi-framework.github.io/Trixi.jl/stable/reference-trixi/#Trixi.DissipationLaxFriedrichsEntropyVariables) 
# for the shallow water system. This is equivalent to `DissipationLocalLaxFriedrichs` if the bottom topography is continuous
@inline function (dissipation::DissipationLaxFriedrichsEntropyVariables)(u_ll, u_rr,
                                                                         orientation_or_normal_direction,
                                                                         equations::ShallowWaterEquations3D)
    h_ll, h_v1_ll, h_v2_ll, h_v3_ll, b_ll = u_ll
    h_rr, h_v1_rr, h_v2_rr, h_v3_rr, b_rr = u_rr
    v1_ll, v2_ll, v3_ll = velocity(u_ll, equations)
    v1_rr, v2_rr, v3_rr = velocity(u_rr, equations)

    # Compute the maximum wave speed
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations)

    # Get some averages
    v1_avg = 0.5 * (v1_ll + v1_rr)
    v2_avg = 0.5 * (v2_ll + v2_rr)
    v3_avg = 0.5 * (v3_ll + v3_rr)

    # Compute dissipation operator
    diss1 = -0.5 * λ * ((h_rr + b_rr) - (h_ll + b_ll))
    diss2 = -0.5 * λ * ((h_v1_rr - h_v1_ll) + v1_avg * (b_rr - b_ll))
    diss3 = -0.5 * λ * ((h_v2_rr - h_v2_ll) + v2_avg * (b_rr - b_ll))
    diss4 = -0.5 * λ * ((h_v3_rr - h_v3_ll) + v3_avg * (b_rr - b_ll))
    return SVector(diss1, diss2, diss3, diss4, zero(eltype(u_ll)))
end

@inline function max_abs_speeds(u, equations::ShallowWaterEquations3D)
    h = waterheight(u, equations)
    v1, v2, v3 = velocity(u, equations)

    c = sqrt(equations.gravity * h)
    return abs(v1) + c, abs(v2) + c, abs(v3) + c
end

# Helper function to extract the velocity vector from the conservative variables
@inline function velocity(u, equations::ShallowWaterEquations3D)
    h, h_v1, h_v2, h_v3, _ = u

    v1 = h_v1 / h
    v2 = h_v2 / h
    v3 = h_v3 / h
    return SVector(v1, v2, v3)
end

# Convert conservative variables to primitive
@inline function cons2prim(u, equations::ShallowWaterEquations3D)
    h, _, _, _, b = u

    H = h + b
    v1, v2, v3 = velocity(u, equations)
    return SVector(H, v1, v2, v3, b)
end

# Convert conservative variables to entropy
# Note, only the first three are the entropy variables, the fourth entry still
# just carries the bottom topography values for convenience
@inline function cons2entropy(u, equations::ShallowWaterEquations3D)
    h, h_v1, h_v2, h_v3, b = u

    v1, v2, v3 = velocity(u, equations)
    v_square = v1^2 + v2^2 + v3^2

    w1 = equations.gravity * (h + b) - 0.5f0 * v_square
    w2 = v1
    w3 = v2
    w4 = v3
    return SVector(w1, w2, w3, w4, b)
end

# Convert entropy variables to conservative
@inline function entropy2cons(w, equations::ShallowWaterEquations3D)
    w1, w2, w3, w4, b = w

    h = (w1 + 0.5f0 * (w2^2 + w3^2 + w4^2)) / equations.gravity - b
    h_v1 = h * w2
    h_v2 = h * w3
    h_v3 = h * w4
    return SVector(h, h_v1, h_v2, h_v3, b)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::ShallowWaterEquations3D)
    H, v1, v2, v3, b = prim

    h = H - b
    h_v1 = h * v1
    h_v2 = h * v2
    h_v3 = h * v3
    return SVector(h, h_v1, h_v2, h_v3, b)
end

@inline function waterheight(u, equations::ShallowWaterEquations3D)
    return u[1]
end

@inline function pressure(u, equations::ShallowWaterEquations3D)
    h = waterheight(u, equations)
    p = 0.5f0 * equations.gravity * h^2
    return p
end

# Entropy function for the shallow water equations is the total energy
@inline function entropy(cons, equations::ShallowWaterEquations3D)
    energy_total(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline function energy_total(cons, equations::ShallowWaterEquations3D)
    h, h_v1, h_v2, h_v3, b = cons

    e = (h_v1^2 + h_v2^2 + h_v3^2) / (2 * h) + 0.5f0 * equations.gravity * h^2 +
        equations.gravity * h * b
    return e
end

# Calculate kinetic energy for a conservative state `cons`
@inline function energy_kinetic(u, equations::ShallowWaterEquations3D)
    h, h_v1, h_v2, h_v3, _ = u
    return (h_v1^2 + h_v2^2 + h_v3^2) / (2 * h)
end

# Calculate potential energy for a conservative state `cons`
@inline function energy_internal(cons, equations::ShallowWaterEquations3D)
    return energy_total(cons, equations) - energy_kinetic(cons, equations)
end

# Convert spherical velocity components to global Cartesian components
@inline function spherical2global(prim, x, ::ShallowWaterEquations3D)
    h, vlat, vlon, vrad, b = prim  # primitive variables

    v1, v2, v3 = spherical2cartesian(vlat, vlon, vrad, x)

    return SVector(h, v1, v2, v3, b)
end
end # @muladd
