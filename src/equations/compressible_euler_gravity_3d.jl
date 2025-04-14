# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent
@doc raw"""
    CompressibleEulerEquationsWithGravity3D(gamma)

The compressible Euler equations with gravity,

References:
- Souza, A. N., He, J., Bischoff, T., Waruszewski, M., Novak, L., Barra, V., ... & Schneider, T. (2023). The flux‐differencing discontinuous galerkin method applied to an idealized fully compressible nonhydrostatic dry atmosphere. Journal of Advances in Modeling Earth Systems, 15(4), e2022MS003527. https://doi.org/10.1029/2022MS003527.
- Waruszewski, M., Kozdon, J. E., Wilcox, L. C., Gibson, T. H., & Giraldo, F. X. (2022). Entropy stable discontinuous Galerkin methods for balance laws in non-conservative form: Applications to the Euler equations with gravity. Journal of Computational Physics, 468, 111507. https://doi.org/10.1016/j.jcp.2022.111507.
"""
struct CompressibleEulerEquationsWithGravity3D{RealT <: Real} <:
       Trixi.AbstractCompressibleEulerEquations{3, 6}
    gamma::RealT               # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications

    function CompressibleEulerEquationsWithGravity3D(gamma)
        γ, inv_gamma_minus_one = promote(gamma, inv(gamma - 1))
        new{typeof(γ)}(γ, inv_gamma_minus_one)
    end
end

Trixi.have_nonconservative_terms(::CompressibleEulerEquationsWithGravity3D) = True()
Trixi.varnames(::typeof(cons2cons), ::CompressibleEulerEquationsWithGravity3D) = ("rho",
                                                                                  "rho_v1",
                                                                                  "rho_v2",
                                                                                  "rho_v3",
                                                                                  "rho_e",
                                                                                  "phi")
Trixi.varnames(::typeof(cons2prim), ::CompressibleEulerEquationsWithGravity3D) = ("rho",
                                                                                  "v1",
                                                                                  "v2",
                                                                                  "v3",
                                                                                  "p",
                                                                                  "phi")

"""
    boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,
                                    equations::CompressibleEulerEquationsWithGravity3D)

Determine the boundary numerical surface flux for a slip wall condition.
Imposes a zero normal velocity at the wall.
Density is taken from the internal solution state and pressure is computed as an
exact solution of a 1D Riemann problem. Further details about this boundary state
are available in the paper:
- J. J. W. van der Vegt and H. van der Ven (2002)
    Slip flow boundary conditions in discontinuous Galerkin discretizations of
    the Euler equations of gas dynamics
    [PDF](https://reports.nlr.nl/bitstream/handle/10921/692/TP-2002-300.pdf?sequence=1)

Details about the 1D pressure Riemann solution can be found in Section 6.3.3 of the book
- Eleuterio F. Toro (2009)
    Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
    3rd edition
    [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)

Should be used together with [`UnstructuredMesh2D`](@ref).
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner,
                                                    normal_direction::AbstractVector,
                                                    x, t,
                                                    surface_flux_function,
                                                    equations::CompressibleEulerEquationsWithGravity3D)
        norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    # Some vector that can't be identical to normal_vector (unless normal_vector == 0)
    tangent1 = SVector(normal_direction[2], normal_direction[3], -normal_direction[1])
    # Orthogonal projection
    tangent1 -= dot(normal, tangent1) * normal
    tangent1 = normalize(tangent1)

    # Third orthogonal vector
    tangent2 = normalize(cross(normal_direction, tangent1))

    # rotate the internal solution state
    u_local = rotate_to_x(u_inner, normal, tangent1, tangent2, equations)

    # compute the primitive variables
    rho_local, v_normal, v_tangent1, v_tangent2, p_local = cons2prim(u_local, equations)

    # Get the solution of the pressure Riemann problem
    # See Section 6.3.3 of
    # Eleuterio F. Toro (2009)
    # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
    # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
    if v_normal <= 0
        sound_speed = sqrt(equations.gamma * p_local / rho_local) # local sound speed
        p_star = p_local *
                 (1 + 0.5f0 * (equations.gamma - 1) * v_normal / sound_speed)^(2 *
                                                                               equations.gamma *
                                                                               equations.inv_gamma_minus_one)
    else # v_normal > 0
        A = 2 / ((equations.gamma + 1) * rho_local)
        B = p_local * (equations.gamma - 1) / (equations.gamma + 1)
        p_star = p_local +
                 0.5f0 * v_normal / A *
                 (v_normal + sqrt(v_normal^2 + 4 * A * (p_local + B)))
    end

    # For the slip wall we directly set the flux as the normal velocity is zero
    return (SVector(zero(eltype(u_inner)),
                    p_star * normal[1],
                    p_star * normal[2],
                    p_star * normal[3],
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner))) * norm_,
            SVector(zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner))) * norm_)
end

# Calculate 3D flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::CompressibleEulerEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, rho_e, phi = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) *
        (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))
    if orientation == 1
        f1 = rho_v1
        f2 = rho_v1 * v1 + p
        f3 = rho_v1 * v2
        f4 = rho_v1 * v3
        f5 = (rho_e + p) * v1
    elseif orientation == 2
        f1 = rho_v2
        f2 = rho_v2 * v1
        f3 = rho_v2 * v2 + p
        f4 = rho_v2 * v3
        f5 = (rho_e + p) * v2
    else
        f1 = rho_v3
        f2 = rho_v3 * v1
        f3 = rho_v3 * v2
        f4 = rho_v3 * v3 + p
        f5 = (rho_e + p) * v3
    end
    return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

# Calculate 2D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::CompressibleEulerEquationsWithGravity3D)
    rho_e = u[5]
    rho, v1, v2, v3, p, _ = cons2prim(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] + v3 * normal_direction[3]
    rho_v_normal = rho * v_normal
    f1 = rho_v_normal
    f2 = rho_v_normal * v1 + p * normal_direction[1]
    f3 = rho_v_normal * v2 + p * normal_direction[2]
    f4 = rho_v_normal * v3 + p * normal_direction[3]
    f5 = (rho_e + p) * v_normal
    return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

"""
    flux_kennedy_gruber(u_ll, u_rr, orientation_or_normal_direction,
                        equations::CompressibleEulerEquationsWithGravity3D)

Kinetic energy preserving two-point flux by
- Kennedy and Gruber (2008)
    Reduced aliasing formulations of the convective terms within the
    Navier-Stokes equations for a compressible fluid
    [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function Trixi.flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::CompressibleEulerEquationsWithGravity3D)
    # Unpack left and right state
    rho_e_ll = u_ll[5]
    rho_e_rr = u_rr[5]
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5 * (rho_ll + rho_rr)
    v1_avg = 0.5 * (v1_ll + v1_rr)
    v2_avg = 0.5 * (v2_ll + v2_rr)
    v3_avg = 0.5 * (v3_ll + v3_rr)
    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2] +
                  v3_avg * normal_direction[3]
    p_avg = 0.5 * (p_ll + p_rr)
    e_avg = 0.5 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * e_avg + p_avg * v_dot_n_avg

    return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

function flux_nonconservative_waruszewski(u_ll, u_rr, normal_direction::AbstractVector,
                                          equations::CompressibleEulerEquationsWithGravity3D)
    rho_ll, _, _, _, _, phi_ll = u_ll
    rho_rr, _, _, _, _, phi_rr = u_rr

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    noncons = ln_mean(rho_ll, rho_rr) * (phi_rr - phi_ll)
    # noncons = 0.5 * (rho_ll + rho_rr) * (phi_rr - phi_ll)

    f0 = zero(eltype(u_ll))
    return SVector(f0, noncons * normal_direction[1], noncons * normal_direction[2],
                   noncons * normal_direction[3], f0, f0)
end

"""
    FluxLMARS(c)(u_ll, u_rr, orientation_or_normal_direction,
                    equations::CompressibleEulerEquationsWithGravity3D)

Low Mach number approximate Riemann solver (LMARS) for atmospheric flows using
an estimate `c` of the speed of sound.

References:
- Xi Chen et al. (2013)
    A Control-Volume Model of the Compressible Euler Equations with a Vertical
    Lagrangian Coordinate
    [DOI: 10.1175/MWR-D-12-00129.1](https://doi.org/10.1175/mwr-d-12-00129.1)
"""
# The struct is already defined in CompressibleEulerEquations2D
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerEquationsWithGravity3D)
    c = flux_lmars.speed_of_sound

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    # Note that this is the same as computing v_ll and v_rr with a normalized normal vector
    # and then multiplying v by `norm_` again, but this version is slightly faster.
    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)
    p = 0.5f0 * (p_ll + p_rr) - 0.5f0 * c * rho * (v_rr - v_ll) / norm_
    v = 0.5f0 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll) * norm_

    # We treat the energy term analogous to the potential temperature term in the paper by
    # Chen et al., i.e. we use p_ll and p_rr, and not p
    if v >= 0
        f1, f2, f3, f4, f5, _ = u_ll * v
        f5 = f5 + p_ll * v
    else
        f1, f2, f3, f4, f5, _ = u_rr * v
        f5 = f5 + p_rr * v
    end

    return SVector(f1,
                   f2 + p * normal_direction[1],
                   f3 + p * normal_direction[2],
                   f4 + p * normal_direction[3],
                   f5, zero(eltype(u_ll)))
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::CompressibleEulerEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1] +
            v2_ll * normal_direction[2] + 
            v3_ll * normal_direction[3])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1] +
            v2_rr * normal_direction[2] +
            v3_rr * normal_direction[3])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
# has been normalized prior to this rotation of the state vector

# Rotate normal vector to x-axis; normal, tangent1 and tangent2 need to be orthonormal
# Called inside `FluxRotated` in `numerical_fluxes.jl` so the directions
# has been normalized prior to this rotation of the state vector
@inline function rotate_to_x(u, normal_vector, tangent1, tangent2,
                             equations::CompressibleEulerEquationsWithGravity3D)
    # Multiply with [ 1   0        0       0   0;
    #                 0   ―  normal_vector ―   0;
    #                 0   ―    tangent1    ―   0;
    #                 0   ―    tangent2    ―   0;
    #                 0   0        0       0   1 ]
    return SVector(u[1],
                   normal_vector[1] * u[2] + normal_vector[2] * u[3] +
                   normal_vector[3] * u[4],
                   tangent1[1] * u[2] + tangent1[2] * u[3] + tangent1[3] * u[4],
                   tangent2[1] * u[2] + tangent2[2] * u[3] + tangent2[3] * u[4],
                   u[5],
                   u[6])
end

@inline function Trixi.max_abs_speeds(u,
                                      equations::CompressibleEulerEquationsWithGravity3D)
    rho, v1, v2, v3, p, _ = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c, abs(v3) + c
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::CompressibleEulerEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, rho_e, phi = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3) - rho * phi)

    return SVector(rho, v1, v2, v3, p, phi)
end

# Convert conservative variables to entropy (see, e.g., Waruszewski et al. (2022))
@inline function Trixi.cons2entropy(u,
                                    equations::CompressibleEulerEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, rho_e, phi = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    v_square = v1^2 + v2^2 + v3^2
    p = (equations.gamma - 1) * (rho_e - 0.5 * rho * v_square - rho * phi)
    s = log(p) - equations.gamma * log(rho)
    rho_p = rho / p

    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         rho_p * (0.5 * v_square - phi)
    w2 = rho_p * v1
    w3 = rho_p * v2
    w4 = rho_p * v3
    w5 = -rho_p

    return SVector(w1, w2, w3, w4, w5, phi)
end

@inline function Trixi.entropy2cons(w,
                                    equations::CompressibleEulerEquationsWithGravity3D)
    # See Waruszewski et al. (2022)
    @unpack gamma = equations

    # convert to entropy `-rho * s`
    # instead of `-rho * s / (gamma - 1)`
    V1, V2, V3, V4, V5, _ = w .* (gamma - 1)
    phi = w[6]

    # s = specific entropy
    V_square = V2^2 + V3^2 + V4^2
    s = gamma - V1 + V_square / (2 * V5) - V5 * phi

    rho_iota = ((gamma - 1) / (-V5)^gamma)^(equations.inv_gamma_minus_one) *
               exp(-s * equations.inv_gamma_minus_one)

    rho = -rho_iota * V5
    rho_v1 = rho_iota * V2
    rho_v2 = rho_iota * V3
    rho_v3 = rho_iota * V4
    rho_e = rho_iota * (1 - V_square / (2 * V5)) + rho * phi
    return SVector(rho, rho_v1, rho_v2, rho_v3, rho_e, phi)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim,
                                 equations::CompressibleEulerEquationsWithGravity3D)
    rho, v1, v2, v3, p, phi = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_v3 = rho * v3
    rho_e = p * equations.inv_gamma_minus_one + 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3) +
            rho * phi
    return SVector(rho, rho_v1, rho_v2, rho_v3, rho_e, phi)
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the
# gravitational potential
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::CompressibleEulerEquationsWithGravity3D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations)
    diss = -0.5 * λ * (u_rr - u_ll)
    return SVector(diss[1], diss[2], diss[3], diss[4], diss[5], zero(eltype(u_ll)))
end

@inline function Trixi.density(u, equations::CompressibleEulerEquationsWithGravity3D)
    rho = u[1]
    return rho
end

@inline function Trixi.velocity(u, equations::CompressibleEulerEquationsWithGravity3D)
    rho = u[1]
    v1 = u[2] / rho
    v2 = u[3] / rho
    v3 = u[4] / rho
    return SVector(v1, v2, v3)
end
end # @muladd
