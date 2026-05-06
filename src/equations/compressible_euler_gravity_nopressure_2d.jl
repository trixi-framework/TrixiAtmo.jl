# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent
@doc raw"""
    CompressibleEulerEquationsWithGravityNoPressure2D(gamma)

The compressible Euler equations with gravity and total energy,
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho v_2 \\ \rho e_{tot} \\ \Phi
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
    \rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ (\rho e_{tot} +p) v_1 \\ 0
\end{pmatrix}
+
\frac{\partial}{\partial y}
\begin{pmatrix}
    \rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ (\rho e_{tot} +p) v_2 \\ 0
\end{pmatrix}
=
\begin{pmatrix}
0 \\ - \rho \nabla \Phi \\ 0 \\ 0
\end{pmatrix}
```
for an ideal gas with ratio of specific heat `gamma` in two space dimensions.

Here, ``\rho`` is the density, ``v_1``,`v_2` the velocities, ``e_{tot}`` the specific total energy, 
which includes the internal, kinetik, and geopotential energy, ``\Phi`` is the gravitational 
geopotential, and
```math
p = (\gamma - 1) \left( \rho e_{tot} - \frac{1}{2} \rho (v_1^2+v_2^2) - \rho \Phi \right)
```
the pressure.

References:
- Souza, A. N., He, J., Bischoff, T., Waruszewski, M., Novak, L., Barra, V., ... & Schneider, T. (2023). The flux‐differencing discontinuous galerkin method applied to an idealized fully compressible nonhydrostatic dry atmosphere. Journal of Advances in Modeling Earth Systems, 15(4), e2022MS003527. https://doi.org/10.1029/2022MS003527.
- Waruszewski, M., Kozdon, J. E., Wilcox, L. C., Gibson, T. H., & Giraldo, F. X. (2022). Entropy stable discontinuous Galerkin methods for balance laws in non-conservative form: Applications to the Euler equations with gravity. Journal of Computational Physics, 468, 111507. https://doi.org/10.1016/j.jcp.2022.111507.
"""
struct CompressibleEulerEquationsWithGravityNoPressure2D{RealT <: Real} <:
       Trixi.AbstractCompressibleEulerEquations{2, 6}
    gamma::RealT               # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications

    function CompressibleEulerEquationsWithGravityNoPressure2D(gamma)
        γ, inv_gamma_minus_one = promote(gamma, inv(gamma - 1))
        new{typeof(γ)}(γ, inv_gamma_minus_one)
    end
end

Trixi.have_nonconservative_terms(::CompressibleEulerEquationsWithGravityNoPressure2D) = True()

# The auxiliary variable is used to store the mean temperature of the element for isothermal-equilibrium-preserving discretizations
Trixi.varnames(::typeof(cons2cons), ::CompressibleEulerEquationsWithGravityNoPressure2D) = ("rho",
                                                                                            "rho_v1",
                                                                                            "rho_v2",
                                                                                            "rho_etot",
                                                                                            "phi",
                                                                                            "aux")
Trixi.varnames(::typeof(cons2prim), ::CompressibleEulerEquationsWithGravityNoPressure2D) = ("rho",
                                                                                            "v1",
                                                                                            "v2",
                                                                                            "p",
                                                                                            "phi",
                                                                                            "aux")

"""
    boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)

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
                                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    # rotate the internal solution state
    u_local = Trixi.rotate_to_x(u_inner, normal, equations)

    # compute the primitive variables
    rho_local, v_normal, v_tangent, p_local, _ = cons2prim(u_local, equations)

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

    # Compute non-conservative term: Since the geopotential is a continuous function, it does not act at the BC
    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    # TODO: compute with surface_flux_function
    # surface_flux, surface_noncons = surface_flux_function
    # surface_noncons(u_ll, u_rr,
    #                                            normal_direction::AbstractVector,
    #                                            equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    noncons = (p_star - p_local)

    # For the slip wall we directly set the flux as the normal velocity is zero
    return (SVector(zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner))) * norm_,
            SVector(zero(eltype(u_inner)),
                    noncons * normal[1],
                    noncons * normal[2],
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner)),
                    zero(eltype(u_inner))) * norm_)
end

"""
    boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                    surface_flux_function, equations::CompressibleEulerEquationsWithGravityNoPressure2D)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation,
                                                    direction, x, t,
                                                    surface_flux_function,
                                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # get the appropriate normal vector from the orientation
    if orientation == 1
        normal_direction = SVector(one(eltype(u_inner)), zero(eltype(u_inner)))
    else # orientation == 2
        normal_direction = SVector(zero(eltype(u_inner)), one(eltype(u_inner)))
    end

    # compute and return the flux using `boundary_condition_slip_wall` routine above
    return boundary_condition_slip_wall(u_inner, normal_direction, direction,
                                        x, t, surface_flux_function, equations)
end

"""
    boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t,
                                    surface_flux_function, equations::CompressibleEulerEquationsWithGravityNoPressure2D)

Should be used together with [`StructuredMesh`](@ref).
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner,
                                                    normal_direction::AbstractVector,
                                                    direction, x, t,
                                                    surface_flux_function,
                                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
    # to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
    if isodd(direction)
        fluxes = boundary_condition_slip_wall(u_inner, -normal_direction,
                                              x, t, surface_flux_function, equations)
        boundary_flux = (-fluxes[1], -fluxes[2])
    else
        boundary_flux = boundary_condition_slip_wall(u_inner, normal_direction,
                                                     x, t, surface_flux_function,
                                                     equations)
    end

    return boundary_flux
end

# Calculate 2D flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, rho_v1, rho_v2, rho_etot, phi, _ = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) *
        (rho_etot - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2) - rho * phi)
    if orientation == 1
        f1 = rho_v1
        f2 = rho_v1 * v1
        f3 = rho_v1 * v2
        f4 = (rho_etot + p) * v1
    else
        f1 = rho_v2
        f2 = rho_v2 * v1
        f3 = rho_v2 * v2
        f4 = (rho_etot + p) * v2
    end
    return SVector(f1, f2, f3, f4, zero(eltype(u)), zero(eltype(u)))
end

# Calculate 2D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_etot = u[4]
    rho, v1, v2, p, _ = cons2prim(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    rho_v_normal = rho * v_normal
    f1 = rho_v_normal
    f2 = rho_v_normal * v1
    f3 = rho_v_normal * v2
    f4 = (rho_etot + p) * v_normal
    return SVector(f1, f2, f3, f4, zero(eltype(u)), zero(eltype(u)))
end

"""
    flux_shima_etal(u_ll, u_rr, orientation_or_normal_direction,
                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)

This flux is is a modification of the original kinetic energy preserving two-point flux by
- Yuichi Kuya, Kosuke Totani and Soshi Kawai (2018)
    Kinetic energy and entropy preserving schemes for compressible flows
    by split convective forms
    [DOI: 10.1016/j.jcp.2018.08.058](https://doi.org/10.1016/j.jcp.2018.08.058)

The modification is in the energy flux to guarantee pressure equilibrium and was developed by
- Nao Shima, Yuichi Kuya, Yoshiharu Tamaki, Soshi Kawai (JCP 2020)
    Preventing spurious pressure oscillations in split convective form discretizations for
    compressible flows
    [DOI: 10.1016/j.jcp.2020.110060](https://doi.org/10.1016/j.jcp.2020.110060)
"""
@inline function Trixi.flux_shima_etal(u_ll, u_rr, orientation::Integer,
                                       equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    kin_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)
    phi_avg = 0.5f0 * (phi_ll + phi_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        pv1_avg = 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
        f1 = rho_avg * v1_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg +
             f1 * phi_avg
    else
        pv2_avg = 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
        f1 = rho_avg * v2_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = p_avg * v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg +
             f1 * phi_avg
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)), zero(eltype(u_ll)))
end

@inline function Trixi.flux_shima_etal(u_ll, u_rr, normal_direction::AbstractVector,
                                       equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr, _ = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)
    phi_avg = 0.5f0 * (phi_ll + phi_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg
    f3 = f1 * v2_avg
    f4 = (f1 * velocity_square_avg +
          p_avg * v_dot_n_avg * equations.inv_gamma_minus_one
          + 0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll)) + f1 * phi_avg

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)), zero(eltype(u_ll)))
end

"""
    flux_kennedy_gruber(u_ll, u_rr, orientation_or_normal_direction,
                        equations::CompressibleEulerEquationsWithGravityNoPressure2D)

Kinetic energy preserving two-point flux by
- Kennedy and Gruber (2008)
    Reduced aliasing formulations of the convective terms within the
    Navier-Stokes equations for a compressible fluid
    [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function Trixi.flux_kennedy_gruber(u_ll, u_rr, orientation::Integer,
                                           equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right state
    rho_etot_ll = u_ll[4]
    rho_etot_rr = u_rr[4]
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    etot_avg = 0.5f0 * (rho_etot_ll / rho_ll + rho_etot_rr / rho_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = rho_avg * v1_avg
        f2 = rho_avg * v1_avg * v1_avg
        f3 = rho_avg * v1_avg * v2_avg
        f4 = (rho_avg * etot_avg + p_avg) * v1_avg
    else
        f1 = rho_avg * v2_avg
        f2 = rho_avg * v2_avg * v1_avg
        f3 = rho_avg * v2_avg * v2_avg
        f4 = (rho_avg * etot_avg + p_avg) * v2_avg
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)), zero(eltype(u_ll)))
end

@inline function Trixi.flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right state
    rho_etot_ll = u_ll[4]
    rho_etot_rr = u_rr[4]
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2]
    p_avg = 0.5f0 * (p_ll + p_rr)
    etot_avg = 0.5f0 * (rho_etot_ll / rho_ll + rho_etot_rr / rho_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg
    f3 = f1 * v2_avg
    f4 = f1 * etot_avg + p_avg * v_dot_n_avg

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)), zero(eltype(u_ll)))
end

"""
    flux_ranocha(u_ll, u_rr, orientation_or_normal_direction,
                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)

Entropy conserving and kinetic energy preserving two-point flux by
- Hendrik Ranocha (2018)
    Generalised Summation-by-Parts Operators and Entropy Stability of Numerical Methods
    for Hyperbolic Balance Laws
    [PhD thesis, TU Braunschweig](https://cuvillier.de/en/shop/publications/7743)
See also
- Hendrik Ranocha (2020)
    Entropy Conserving and Kinetic Energy Preserving Numerical Methods for
    the Euler Equations Using Summation-by-Parts Operators
    [Proceedings of ICOSAHOM 2018](https://doi.org/10.1007/978-3-030-39647-3_42)
"""
@inline function Trixi.flux_ranocha(u_ll, u_rr, orientation::Integer,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
    # in exact arithmetic since
    #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
    #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)
    phi_avg = 0.5f0 * (phi_ll + phi_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = rho_mean * v1_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = f1 *
             (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
             0.5f0 * (p_ll * v1_rr + p_rr * v1_ll) + f1 * phi_avg
    else
        f1 = rho_mean * v2_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = f1 *
             (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
             0.5f0 * (p_ll * v2_rr + p_rr * v2_ll) + f1 * phi_avg
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)), zero(eltype(u_ll)))
end

@inline function Trixi.flux_ranocha(u_ll, u_rr, normal_direction::AbstractVector,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr, _ = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
    # in exact arithmetic since
    #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
    #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)
    phi_avg = 0.5f0 * (phi_ll + phi_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg
    f3 = f1 * v2_avg
    f4 = (f1 * (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one)
          +
          0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll)
          + f1 * phi_avg)

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)), zero(eltype(u_ll)))
end

function flux_nonconservative_waruszewski_etal(u_ll, u_rr,
                                               normal_direction::AbstractVector,
                                               equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, _, _, p_ll, phi_ll, _ = cons2prim(u_ll, equations)
    rho_rr, _, _, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    p_jump = (p_rr - p_ll)
    noncons = p_jump + ln_mean(rho_ll, rho_rr) * (phi_rr - phi_ll)

    f0 = zero(eltype(u_ll))
    return SVector(f0, noncons * normal_direction[1], noncons * normal_direction[2],
                   f0, f0, f0)
end

function flux_nonconservative_waruszewski_etal(u_ll, u_rr, orientation::Integer,
                                               equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, _, _, p_ll, phi_ll, _ = cons2prim(u_ll, equations)
    rho_rr, _, _, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    p_jump = (p_rr - p_ll)
    noncons = p_jump + ln_mean(rho_ll, rho_rr) * (phi_rr - phi_ll)

    f0 = zero(eltype(u_ll))
    if orientation == 1
        return SVector(f0, noncons, f0, f0, f0, f0)
    else #if orientation == 2
        return SVector(f0, f0, noncons, f0, f0, f0)
    end
end

# For `VolumeIntegralSubcellLimiting` the nonconservative flux is created as a callable struct to 
# enable dispatch on the type of the nonconservative term (local / jump).
struct FluxNonConservativeChandrashekarIsothermal <:
       Trixi.FluxNonConservative{Trixi.NonConservativeJump()}
end

Trixi.n_nonconservative_terms(::FluxNonConservativeChandrashekarIsothermal) = 2

const flux_nonconservative_chandrashekar_isothermal = FluxNonConservativeChandrashekarIsothermal()

@inline function (noncons_flux::FluxNonConservativeChandrashekarIsothermal)(u_ll, u_rr,
                                                                            normal_direction::AbstractVector,
                                                                            equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, _, _, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)
    _, _, _, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # We use the mean element temperature on the left to keep the iso-thermal correction
    # completely element-local, as in:
    # - Chandrashekar, P., & Zenk, M. (2017). Well-balanced nodal discontinuous Galerkin method for Euler 
    #   equations with gravity. Journal of Scientific Computing, 71(3), 1062-1093.
    RT = RT_ll

    e_ll = exp(phi_ll / RT)
    e_rr = exp(phi_rr / RT)

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    noncons = -rho_ll * RT * e_ll * (1 / e_rr - 1 / e_ll) + (p_rr - p_ll)

    f0 = zero(eltype(u_ll))
    return SVector(f0, noncons * normal_direction[1], noncons * normal_direction[2], f0,
                   f0, f0)
end

@inline function (noncons_flux::FluxNonConservativeChandrashekarIsothermal)(u_ll, u_rr,
                                                                            orientation::Integer,
                                                                            equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, _, _, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)
    _, _, _, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # We use the mean element temperature on the left to keep the iso-thermal correction
    # completely element-local, as in:
    # - Chandrashekar, P., & Zenk, M. (2017). Well-balanced nodal discontinuous Galerkin method for Euler 
    #   equations with gravity. Journal of Scientific Computing, 71(3), 1062-1093.
    RT = RT_ll

    e_ll = exp(phi_ll / RT)
    e_rr = exp(phi_rr / RT)

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    noncons = -rho_ll * RT * e_ll * (1 / e_rr - 1 / e_ll) + (p_rr - p_ll)

    f0 = zero(eltype(u_ll))
    if orientation == 1
        return SVector(f0, noncons, f0, f0, f0, f0)
    else #if orientation == 2
        return SVector(f0, f0, noncons, f0, f0, f0)
    end
end

@inline function flux_nonconservative_chandrashekar_isothermal(u_ll,
                                                               orientation::Integer,
                                                               equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                               nonconservative_type::Trixi.NonConservativeLocal,
                                                               nonconservative_term::Integer)

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    if nonconservative_term == 1
        rho_ll, v1_ll, v2_ll, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)
        e_ll = exp(phi_ll / RT_ll)

        flux = -rho_ll * RT_ll * e_ll

        f0 = zero(eltype(u_ll))
        if orientation == 1
            return SVector(f0, flux, f0, f0, f0, f0)
        else #if orientation == 2
            return SVector(f0, f0, flux, f0, f0, f0)
        end
    else #if nonconservative_term == 2
        flux = 1

        f0 = zero(eltype(u_ll))
        if orientation == 1
            return SVector(f0, flux, f0, f0, f0, f0)
        else #if orientation == 2
            return SVector(f0, f0, flux, f0, f0, f0)
        end
    end
end

@inline function flux_nonconservative_chandrashekar_isothermal(u_ll,
                                                               normal_direction::AbstractVector,
                                                               equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                               nonconservative_type::Trixi.NonConservativeLocal,
                                                               nonconservative_term::Integer)
    rho_ll, _, _, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)

    f0 = zero(eltype(u_ll))

    if nonconservative_term == 1
        # We use the mean element temperature on the left to keep the iso-thermal correction
        # completely element-local, as in:
        # - Chandrashekar, P., & Zenk, M. (2017). Well-balanced nodal discontinuous Galerkin method for Euler 
        #   equations with gravity. Journal of Scientific Computing, 71(3), 1062-1093.
        RT = RT_ll

        e_ll = exp(phi_ll / RT)

        # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
        noncons = -rho_ll * RT * e_ll
    else # nonconservative_term == 2
        noncons = one(eltype(u_ll))
    end

    return SVector(f0, noncons * normal_direction[1], noncons * normal_direction[2], f0,
                   f0, f0)
end

@inline function flux_nonconservative_chandrashekar_isothermal(u_ll, u_rr,
                                                               orientation::Integer,
                                                               equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                               nonconservative_type::Trixi.NonConservativeJump,
                                                               nonconservative_term::Integer)
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    # We use the mean element temperature on the left to keep the iso-thermal correction
    # completely element-local, as in:
    # - Chandrashekar, P., & Zenk, M. (2017). Well-balanced nodal discontinuous Galerkin method for Euler 
    #   equations with gravity. Journal of Scientific Computing, 71(3), 1062-1093.
    RT = RT_ll

    # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
    if nonconservative_term == 1
        inv_e_ll = exp(-phi_ll / RT)
        inv_e_rr = exp(-phi_rr / RT)

        flux = (inv_e_rr - inv_e_ll)

        f0 = zero(eltype(u_ll))
        if orientation == 1
            return SVector(f0, flux, f0, f0, f0, f0)
        else #if orientation == 2
            return SVector(f0, f0, flux, f0, f0, f0)
        end
    else #if nonconservative_term == 2
        flux = p_rr - p_ll

        f0 = zero(eltype(u_ll))
        if orientation == 1
            return SVector(f0, flux, f0, f0, f0, f0)
        else #if orientation == 2
            return SVector(f0, f0, flux, f0, f0, f0)
        end
    end
end

@inline function flux_nonconservative_chandrashekar_isothermal(u_ll, u_rr,
                                                               normal_direction::AbstractVector,
                                                               equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                               nonconservative_type::Trixi.NonConservativeJump,
                                                               nonconservative_term::Integer)
    rho_ll, _, _, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)
    _, _, _, p_rr, phi_rr, _ = cons2prim(u_rr, equations)

    f0 = zero(eltype(u_ll))

    if nonconservative_term == 1
        # We use the mean element temperature on the left to keep the iso-thermal correction
        # completely element-local, as in:
        # - Chandrashekar, P., & Zenk, M. (2017). Well-balanced nodal discontinuous Galerkin method for Euler 
        #   equations with gravity. Journal of Scientific Computing, 71(3), 1062-1093.
        RT = RT_ll

        e_ll = exp(phi_ll / RT)
        e_rr = exp(phi_rr / RT)

        # We omit the 0.5 in the density average since Trixi.jl always multiplies the non-conservative flux with 0.5
        noncons = (1 / e_rr - 1 / e_ll)

    else # nonconservative_term == 2
        noncons = (p_rr - p_ll)
    end

    return SVector(f0, noncons, noncons, f0, f0, f0)
end

"""
    FluxLMARS(c)(u_ll, u_rr, orientation_or_normal_direction,
                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)

Low Mach number approximate Riemann solver (LMARS) for atmospheric flows using
an estimate `c` of the speed of sound.

References:
- Xi Chen et al. (2013)
    A Control-Volume Model of the Compressible Euler Equations with a Vertical
    Lagrangian Coordinate
    [DOI: 10.1175/MWR-D-12-00129.1](https://doi.org/10.1175/mwr-d-12-00129.1)
"""
# The struct is already defined in CompressibleEulerEquations2D

# We add the "upwinding" pressure terms here, but not the average pressure terms as in CompressibleEulerEquationsWithGravity2D
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
                                         equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    c = flux_lmars.speed_of_sound

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end

    rho = 0.5f0 * (rho_ll + rho_rr)
    p = -0.5f0 * c * rho * (v_rr - v_ll)
    v = 0.5f0 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll)

    # We treat the energy term analogous to the potential temperature term in the paper by
    # Chen et al., i.e. we use p_ll and p_rr, and not p
    if v >= 0
        f1, f2, f3, f4, _ = v * u_ll
        f4 = f4 + p_ll * v
    else
        f1, f2, f3, f4, _ = v * u_rr
        f4 = f4 + p_rr * v
    end

    if orientation == 1
        f2 = f2 + p
    else # orientation == 2
        f3 = f3 + p
    end
    f0 = zero(eltype(u_ll))
    return SVector(f1, f2, f3, f4, f0, f0)
end

# We add the "upwinding" pressure terms here, but not the average pressure terms as in CompressibleEulerEquationsWithGravity2D
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    c = flux_lmars.speed_of_sound

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Note that this is the same as computing v_ll and v_rr with a normalized normal vector
    # and then multiplying v by `norm_` again, but this version is slightly faster.
    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)
    p = -0.5f0 * c * rho * (v_rr - v_ll) / norm_
    v = 0.5f0 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll) * norm_

    # We treat the energy term analogous to the potential temperature term in the paper by
    # Chen et al., i.e. we use p_ll and p_rr, and not p
    if v >= 0
        f1, f2, f3, f4, _ = u_ll * v
        f4 = f4 + p_ll * v
    else
        f1, f2, f3, f4, _ = u_rr * v
        f4 = f4 + p_rr * v
    end
    f0 = zero(eltype(u_ll))
    return SVector(f1,
                   f2 + p * normal_direction[1],
                   f3 + p * normal_direction[2],
                   f4, f0, f0)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    λ_max = max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1]
            +
            v2_ll * normal_direction[2])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            +
            v2_rr * normal_direction[2])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Calculate minimum and maximum wave speeds for HLL-type fluxes
@inline function Trixi.min_max_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    if orientation == 1 # x-direction
        λ_min = v1_ll - sqrt(equations.gamma * p_ll / rho_ll)
        λ_max = v1_rr + sqrt(equations.gamma * p_rr / rho_rr)
    else # y-direction
        λ_min = v2_ll - sqrt(equations.gamma * p_ll / rho_ll)
        λ_max = v2_rr + sqrt(equations.gamma * p_rr / rho_rr)
    end

    return λ_min, λ_max
end

@inline function Trixi.min_max_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

    v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    norm_ = norm(normal_direction)
    # The v_normals are already scaled by the norm
    λ_min = v_normal_ll - sqrt(equations.gamma * p_ll / rho_ll) * norm_
    λ_max = v_normal_rr + sqrt(equations.gamma * p_rr / rho_rr) * norm_

    return λ_min, λ_max
end

# Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
# has been normalized prior to this rotation of the state vector
@inline function Trixi.rotate_to_x(u, normal_vector,
                                   equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    # Apply the 2D rotation matrix with normal and tangent directions of the form
    # [ 1    0    0   0;
    #   0   n_1  n_2  0;
    #   0   t_1  t_2  0;
    #   0    0    0   1 ]
    # where t_1 = -n_2 and t_2 = n_1

    return SVector(u[1],
                   c * u[2] + s * u[3],
                   -s * u[2] + c * u[3],
                   u[4],
                   u[5],
                   u[6])
end

# Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
# has been normalized prior to this back-rotation of the state vector
@inline function Trixi.rotate_from_x(u, normal_vector,
                                     equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    # Apply the 2D back-rotation matrix with normal and tangent directions of the form
    # [ 1    0    0   0;
    #   0   n_1  t_1  0;
    #   0   n_2  t_2  0;
    #   0    0    0   1 ]
    # where t_1 = -n_2 and t_2 = n_1

    return SVector(u[1],
                   c * u[2] - s * u[3],
                   s * u[2] + c * u[3],
                   u[4],
                   u[5],
                   u[6])
end

@inline function Trixi.max_abs_speeds(u,
                                      equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, v1, v2, p, _ = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u,
                                 equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, rho_v1, rho_v2, rho_etot, phi, aux = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) *
        (rho_etot - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2) - rho * phi)

    return SVector(rho, v1, v2, p, phi, aux)
end

# Convert conservative variables to entropy (see, e.g., Waruszewski et al. (2022))
@inline function Trixi.cons2entropy(u,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, rho_v1, rho_v2, rho_etot, phi, aux = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v_square = v1^2 + v2^2
    p = (equations.gamma - 1) * (rho_etot - 0.5f0 * rho * v_square - rho * phi)
    s = log(p) - equations.gamma * log(rho)
    rho_p = rho / p

    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         rho_p * (0.5f0 * v_square - phi)
    w2 = rho_p * v1
    w3 = rho_p * v2
    w4 = -rho_p

    return SVector(w1, w2, w3, w4, phi, aux)
end

@inline function Trixi.entropy2cons(w,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # See Waruszewski et al. (2022)
    @unpack gamma = equations

    # convert to entropy `-rho * s`
    # instead of `-rho * s / (gamma - 1)`
    V1, V2, V3, V5, _ = w .* (gamma - 1)
    phi = w[5]
    aux = w[6]

    # s = specific entropy
    s = gamma - V1 + (V2^2 + V3^2) / (2 * V5) - V5 * phi

    rho_iota = ((gamma - 1) / (-V5)^gamma)^(equations.inv_gamma_minus_one) *
               exp(-s * equations.inv_gamma_minus_one)

    rho = -rho_iota * V5
    rho_v1 = rho_iota * V2
    rho_v2 = rho_iota * V3
    rho_etot = rho_iota * (1 - (V2^2 + V3^2) / (2 * V5)) + rho * phi
    return SVector(rho, rho_v1, rho_v2, rho_etot, phi, aux)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim,
                                 equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, v1, v2, p, phi, aux = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_etot = p * equations.inv_gamma_minus_one + 0.5f0 * (rho_v1 * v1 + rho_v2 * v2) +
               rho * phi
    return SVector(rho, rho_v1, rho_v2, rho_etot, phi, aux)
end

@inline function Trixi.density(u,
                               equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho = u[1]
    return rho
end

@inline function Trixi.pressure(u,
                                equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, rho_v1, rho_v2, rho_etot, phi, _ = u
    p = (equations.gamma - 1) *
        (rho_etot - 0.5f0 * (rho_v1^2 + rho_v2^2) / rho - rho * phi)
    return p
end

@inline function Trixi.density_pressure(u,
                                        equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, rho_v1, rho_v2, rho_etot, phi, _ = u
    rho_times_p = (equations.gamma - 1) *
                  (rho * rho_etot - 0.5f0 * (rho_v1^2 + rho_v2^2) - rho^2 * phi)
    return rho_times_p
end

# Calculate thermodynamic entropy for a conservative state `cons`
@inline function entropy_thermodynamic(cons,
                                       equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Pressure
    p = (equations.gamma - 1) *
        (cons[4] - 0.5f0 * (cons[2]^2 + cons[3]^2) / cons[1] - cons[1] * cons[5])

    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])

    return s
end

# Calculate mathematical entropy for a conservative state `cons`
@inline function entropy_math(cons,
                              equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Mathematical entropy
    S = -entropy_thermodynamic(cons, equations) * cons[1] *
        equations.inv_gamma_minus_one

    return S
end

# Default entropy is the mathematical entropy
@inline Trixi.entropy(cons, equations::CompressibleEulerEquationsWithGravityNoPressure2D) = entropy_math(cons,
                                                                                                         equations)

# Calculate total energy for a conservative state `cons`
@inline Trixi.energy_total(cons, ::CompressibleEulerEquationsWithGravityNoPressure2D) = cons[4]

# Calculate kinetic energy for a conservative state `cons`
@inline function Trixi.energy_kinetic(u,
                                      equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho, rho_v1, rho_v2, rho_etot, _ = u
    return (rho_v1^2 + rho_v2^2) / (2 * rho)
end

# Calculate internal energy for a conservative state `cons`
@inline function Trixi.energy_internal(cons,
                                       equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    return energy_total(cons, equations) - energy_kinetic(cons, equations) -
           cons[1] * cons[5]
end

@inline function Trixi.velocity(u,
                                equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho = u[1]
    v1 = u[2] / rho
    v2 = u[3] / rho
    return SVector(v1, v2)
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the
# gravitational potential
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations)
    diss = -0.5f0 * λ * (u_rr - u_ll)
    f0 = zero(eltype(u_ll))
    return SVector(diss[1], diss[2], diss[3], diss[4], f0, f0)
end

# State validation for Newton-bisection method of subcell IDP limiting
@inline function Base.isvalid(u,
                              equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    if u[1] <= 0 || pressure(u, equations) <= 0
        return false
    end
    return true
end

#######################################################################
# Hydrostatic reconstruction

struct FluxHydrostaticReconstruction{NumericalFlux, HydrostaticReconstruction}
    numerical_flux::NumericalFlux
    hydrostatic_reconstruction::HydrostaticReconstruction
end

@inline function (numflux::FluxHydrostaticReconstruction)(u_ll, u_rr,
                                                          orientation_or_normal_direction,
                                                          equations::Trixi.AbstractEquations)
    @unpack numerical_flux, hydrostatic_reconstruction = numflux

    # Create the reconstructed left/right solution states in conservative form
    u_ll_star, u_rr_star = hydrostatic_reconstruction(u_ll, u_rr, equations)

    # Use the reconstructed states to compute the numerical surface flux
    return numerical_flux(u_ll_star, u_rr_star, orientation_or_normal_direction,
                          equations)
end

# Hydrostatic reconstruction for the isothermal equilibrium from the paper:
# Ziming Chen, Yingjuan Zhang, Gang Li, Shouguo Qian (2022)
# "A well-balanced Runge-Kutta discontinuous Galerkin method for the Euler equations in isothermal 
# hydrostatic state under gravitational field"
# [DOI:10.1016/j.camwa.2022.05.025](https://doi.org/10.1016/j.camwa.2022.05.025)
@inline function hydrostatic_reconstruction_isothermal(u_ll, u_rr,
                                            equations::Union{CompressibleEulerEquationsWithGravityNoPressure2D,
                                                             CompressibleEulerEquationsWithGravity2D})
    # Unpack left and right states
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll, RT_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr, RT_rr = cons2prim(u_rr, equations)

    rho_e_ll = u_ll[4]
    rho_e_rr = u_rr[4]

    psi_ll = p_ll / rho_ll * log(rho_ll)
    psi_rr = p_rr / rho_rr * log(rho_rr)

    # Compute equilibrium potential (# TODO: For general case we need to add phi_num - phi_exact))
    psi_eq_ll = psi_ll # phi_ll - phi_exact_ll
    psi_eq_rr = psi_rr # phi_rr - phi_exact_rr

    # We use the mean element temperature on the left to keep the iso-thermal correction
    # completely element-local, as in:
    # - Chandrashekar, P., & Zenk, M. (2017). Well-balanced nodal discontinuous Galerkin method for Euler 
    #   equations with gravity. Journal of Scientific Computing, 71(3), 1062-1093.
    RT0 = RT_ll

    # Compute equlibrium state
    rho_eq_ll = exp(psi_eq_ll / (RT0))
    rho_eq_rr = exp(psi_eq_rr / (RT0))
    p_eq_ll = RT0 * rho_eq_ll
    p_eq_rr = RT0 * rho_eq_rr
    rho_e_eq_ll = p_eq_ll / (equations.gamma - 1) + rho_eq_ll * phi_ll
    rho_e_eq_rr = p_eq_rr / (equations.gamma - 1) + rho_eq_rr * phi_rr

    u_eq_ll = SVector(rho_eq_ll, 0.0, 0.0, rho_e_eq_ll, phi_ll, RT_ll)
    u_eq_rr = SVector(rho_eq_rr, 0.0, 0.0, rho_e_eq_rr, phi_rr, RT_rr)

    # Compute residual contribution
    u_res_ll = u_ll - u_eq_ll
    u_res_rr = u_rr - u_eq_rr

    # Reconstruct phi
    phi_star = max(phi_ll, phi_rr)

    # Compute reconstructed equilibrium potential
    psi_eq_ll = psi_ll + phi_ll - phi_star
    psi_eq_rr = psi_rr + phi_rr - phi_star

    # Compute reconstructed equilibrium state
    rho_eq_ll = exp(psi_eq_ll / (RT0))
    rho_eq_rr = exp(psi_eq_rr / (RT0))
    p_eq_ll = RT0 * rho_eq_ll
    p_eq_rr = RT0 * rho_eq_rr
    rho_e_eq_ll = p_eq_ll / (equations.gamma - 1) + rho_eq_ll * phi_star
    rho_e_eq_rr = p_eq_rr / (equations.gamma - 1) + rho_eq_rr * phi_star

    u_eq_ll = SVector(rho_eq_ll, 0.0, 0.0, rho_e_eq_ll, phi_star, RT_ll)
    u_eq_rr = SVector(rho_eq_rr, 0.0, 0.0, rho_e_eq_rr, phi_star, RT_rr)

    # Compute reconstructed state
    u_star_ll = u_eq_ll + u_res_ll
    u_star_rr = u_eq_rr + u_res_rr

    return u_star_ll, u_star_rr
end

# Get outer state for the min/max limiters
@inline function Trixi.get_boundary_outer_state(u_inner, t,
                                                boundary_condition::typeof(boundary_condition_slip_wall),
                                                orientation, boundary_index,
                                                mesh::TreeMesh{2},
                                                equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                dg, cache, indices...)
    if orientation == 1
        #if boundary_index == 1 # Element is on the right, boundary on the left
        return SVector(u_inner[1],
                       u_inner[2] - 2 * u_inner[2],
                       u_inner[3],
                       u_inner[4],
                       u_inner[5],
                       u_inner[6])
    else
        return SVector(u_inner[1],
                       u_inner[2],
                       u_inner[3] - 2 * u_inner[3],
                       u_inner[4],
                       u_inner[5],
                       u_inner[6])
    end
end

@inline function Trixi.get_boundary_outer_state(u_inner, t,
                                                boundary_condition::typeof(boundary_condition_slip_wall),
                                                normal_direction::AbstractVector,
                                                mesh::P4estMesh{2},
                                                equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                dg, cache, indices...)
    factor = (normal_direction[1] * u_inner[2] + normal_direction[2] * u_inner[3])
    u_normal = (factor / sum(normal_direction .^ 2)) * normal_direction

    return SVector(u_inner[1],
                   u_inner[2] - 2 * u_normal[1],
                   u_inner[3] - 2 * u_normal[2],
                   u_inner[4],
                   u_inner[5],
                   u_inner[6])
end
end # @muladd
