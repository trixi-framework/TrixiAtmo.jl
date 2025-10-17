# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    CompressibleEulerInternalKineticEnergyEquationsWithGravity2D(gamma)

The compressible Euler equations with gravity
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho v_2 \\ \rho e
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
 \rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ (\rho e +p) v_1
\end{pmatrix}
+
\frac{\partial}{\partial y}
\begin{pmatrix}
\rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ (\rho e +p) v_2
\end{pmatrix}
=
\begin{pmatrix}
0 \\ - rho \frac{\partial}{\partial x} \phi \\ - rho \frac{\partial}{\partial y} \phi \\ - rho v \frac{\partial}{\partial x} \phi - rho v \frac{\partial}{\partial y} \phi
\end{pmatrix}
```
for an ideal gas with ratio of specific heats `gamma`
in two space dimensions.
Here, ``\rho`` is the density, ``v_1``, ``v_2`` the velocities, ``e`` the specific total energy **rather than** specific internal energy, ``phi`` is the potential and
```math
p = (\gamma - 1) \left( \rho e - \frac{1}{2} \rho (v_1^2+v_2^2) \right)
```
the pressure.
"""
struct CompressibleEulerInternalKineticEnergyEquationsWithGravity2D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{2, 5}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats 
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    function CompressibleEulerInternalKineticEnergyEquationsWithGravity2D(; c_p, c_v,
                                                                          gravity)
        c_p, c_v, g = promote(c_p, c_v, gravity)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        return new{typeof(c_p)}(p_0, c_p, c_v, g, R,
                                gamma,
                                inv_gamma_minus_one)
    end
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    ("rho", "rho_v1", "rho_v2", "rho_e", "phi")
end
varnames(::typeof(cons2prim), ::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D) = ("rho",
                                                                                                 "v1",
                                                                                                 "v2",
                                                                                                 "p",
                                                                                                 "phi")

have_nonconservative_terms(::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D) = Trixi.True()

# Slip-wall boundary condition
# Determine the boundary numerical surface flux for a slip wall condition.
# Imposes a zero normal velocity at the wall.
@inline function boundary_condition_slip_wall(u_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_functions,
                                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4], u_inner[5])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

@inline function boundary_condition_slip_wall(u_inner, orientation,
                                              direction, x, t,
                                              surface_flux_functions,
                                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4],
                             u_inner[5])
    else # orientation == 2
        u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4],
                             u_inner[5])
    end

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
                                                     equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
                                                     equations)
    end

    return flux, noncons_flux
end

# Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer,
                      equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2))
    if orientation == 1
        f1 = rho_v1
        f2 = rho_v1 * v1 + p
        f3 = rho_v1 * v2
        f4 = (rho_e + p) * v1
    else
        f1 = rho_v2
        f2 = rho_v2 * v1
        f3 = rho_v2 * v2 + p
        f4 = (rho_e + p) * v2
    end
    return SVector(f1, f2, f3, f4, zero(eltype(u)))
end

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_e = last(u)
    rho, v1, v2, p = cons2prim(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    rho_v_normal = rho * v_normal
    f1 = rho_v_normal
    f2 = rho_v_normal * v1 + p * normal_direction[1]
    f3 = rho_v_normal * v2 + p * normal_direction[2]
    f4 = (rho_e + p) * v_normal
    return SVector(f1, f2, f3, f4, zero(eltype(u)))
end

"""
	flux_nonconservative_waruzewski_etal(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)

Well-balanced gravity term for isothermal background state
-  Maciej Waruszewski and Jeremy E. Kozdon and Lucas C. Wilcox and Thomas H. Gibson and Francis X. Giraldo (2022)
   Entropy stable discontinuous {G}alerkin methods for balance laws 
   in non-conservative form: Applications to the {E}uler equations with gravity
   [DOI: 10.1016/j.jcp.2022.111507](https://doi.org/10.1016/j.jcp.2022.111507)

The well balanced on curvilinear coordinates was proven by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_waruzewski_etal(u_ll, u_rr,
                                                      normal_direction::AbstractVector,
                                                      equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, rho_v1_ll, rho_v2_ll, _, phi_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, _, phi_rr = u_rr
    v1_ll = rho_v1_ll / rho_ll
    v1_rr = rho_v1_rr / rho_rr
    v2_ll = rho_v2_ll / rho_ll
    v2_rr = rho_v2_rr / rho_rr
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    rho_avg = ln_mean(rho_ll, rho_rr)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   normal_direction[1] * rho_avg * v1_avg * phi_jump +
                   normal_direction[2] * rho_avg * v2_avg * phi_jump,
                   zero(eltype(u_ll)))
end

"""
	flux_nonconservative_artiano_etal(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)

Well-balanced gravity term for constant potential temperature background state by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations 
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_artiano_etal(u_ll, u_rr,
                                                   normal_direction::AbstractVector,
                                                   equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, rho_v1_ll, rho_v2_ll, _, phi_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, _, phi_rr = u_rr
    v1_ll = rho_v1_ll / rho_ll
    v1_rr = rho_v1_rr / rho_rr
    v2_ll = rho_v2_ll / rho_ll
    v2_rr = rho_v2_rr / rho_rr
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   normal_direction[1] * rho_avg * v1_avg * phi_jump +
                   normal_direction[2] * rho_avg * v2_avg * phi_jump,
                   zero(eltype(u_ll)))
end

"""
	flux_nonconservative_souza_etal(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)

-  Souza et al. 
   The Flux-Differencing Discontinuous {G}alerkin Method Applied to 
   an Idealized Fully Compressible Nonhydrostatic Dry Atmosphere
   [DOI: 10.1029/2022MS003527] (https://doi.org/10.1029/2022MS003527)
"""
@inline function flux_nonconservative_souza_etal(u_ll, u_rr,
                                                 normal_direction::AbstractVector,
                                                 equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, rho_v1_ll, rho_v2_ll, _, phi_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, _, phi_rr = u_rr
    v1_ll = rho_v1_ll / rho_ll
    v1_rr = rho_v1_rr / rho_rr
    v2_ll = rho_v2_ll / rho_ll
    v2_rr = rho_v2_rr / rho_rr
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   normal_direction[1] * rho_avg * v1_avg * phi_jump +
                   normal_direction[2] * rho_avg * v2_avg * phi_jump,
                   zero(eltype(u_ll)))
end

"""
    flux_shima_etal(u_ll, u_rr, orientation_or_normal_direction,
                    equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)

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
@inline function flux_shima_etal(u_ll, u_rr, orientation::Integer,
                                 equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    kin_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        pv1_avg = 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
        f1 = rho_avg * v1_avg
        f2 = f1 * v1_avg + p_avg
        f3 = f1 * v2_avg
        f4 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg
    else
        pv2_avg = 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
        f1 = rho_avg * v2_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg + p_avg
        f4 = p_avg * v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function flux_shima_etal(u_ll, u_rr, normal_direction::AbstractVector,
                                 equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = (f1 * velocity_square_avg +
          p_avg * v_dot_n_avg * equations.inv_gamma_minus_one
          + 0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll))

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

"""
    flux_kennedy_gruber(u_ll, u_rr, orientation_or_normal_direction,
                        equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)

Kinetic energy preserving two-point flux by
- Kennedy and Gruber (2008)
  Reduced aliasing formulations of the convective terms within the
  Navier-Stokes equations for a compressible fluid
  [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function flux_kennedy_gruber(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_e_ll = last(u_ll)
    rho_e_rr = last(u_rr)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = rho_avg * v1_avg
        f2 = rho_avg * v1_avg * v1_avg + p_avg
        f3 = rho_avg * v1_avg * v2_avg
        f4 = (rho_avg * e_avg + p_avg) * v1_avg
    else
        f1 = rho_avg * v2_avg
        f2 = rho_avg * v2_avg * v1_avg
        f3 = rho_avg * v2_avg * v2_avg + p_avg
        f4 = (rho_avg * e_avg + p_avg) * v2_avg
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_e_ll = last(u_ll)
    rho_e_rr = last(u_rr)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2]
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * e_avg + p_avg * v_dot_n_avg

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

"""
    flux_ranocha(u_ll, u_rr, orientation_or_normal_direction,
                 equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)

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
@inline function flux_ranocha(u_ll, u_rr, orientation::Integer,
                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
    # in exact arithmetic since
    #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
    #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = rho_mean * v1_avg
        f2 = f1 * v1_avg + p_avg
        f3 = f1 * v2_avg
        f4 = f1 *
             (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
             0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
    else
        f1 = rho_mean * v2_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg + p_avg
        f4 = f1 *
             (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
             0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function flux_ranocha(u_ll, u_rr, normal_direction::AbstractVector,
                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
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
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = (f1 * (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one)
          +
          0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll))

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
                                         equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    c = flux_lmars.speed_of_sound

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end

    rho = 0.5f0 * (rho_ll + rho_rr)
    p = 0.5f0 * (p_ll + p_rr) - 0.5f0 * c * rho * (v_rr - v_ll)
    v = 0.5f0 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll)

    # We treat the energy term analogous to the potential temperature term in the paper by
    # Chen et al., i.e. we use p_ll and p_rr, and not p
    if v >= 0
        f1, f2, f3, f4 = v * u_ll
        f4 = f4 + p_ll * v
    else
        f1, f2, f3, f4 = v * u_rr
        f4 = f4 + p_rr * v
    end

    if orientation == 1
        f2 = f2 + p
    else # orientation == 2
        f3 = f3 + p
    end

    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    c = flux_lmars.speed_of_sound

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Note that this is the same as computing v_ll and v_rr with a normalized normal vector
    # and then multiplying v by `norm_` again, but this version is slightly faster.
    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)
    p = 0.5f0 * (p_ll + p_rr) - 0.5f0 * c * rho * (v_rr - v_ll) / norm_
    v = 0.5f0 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll) * norm_

    # We treat the energy term analogous to the potential temperature term in the paper by
    # Chen et al., i.e. we use p_ll and p_rr, and not p
    if v >= 0
        f1, f2, f3, f4 = u_ll * v
        f4 = f4 + p_ll * v
    else
        f1, f2, f3, f4 = u_rr * v
        f4 = f4 + p_rr * v
    end

    return SVector(f1,
                   f2 + p * normal_direction[1],
                   f3 + p * normal_direction[2],
                   f4, zero(eltype(u_ll)))
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

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

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

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

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, orientation::Integer,
                               equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

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

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speeds
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

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_,
               abs(v_rr) + c_rr * norm_)
end

# Calculate estimate for minimum and maximum wave speeds for HLL-type fluxes
@inline function min_max_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    if orientation == 1 # x-direction
        λ_min = v1_ll - sqrt(equations.gamma * p_ll / rho_ll)
        λ_max = v1_rr + sqrt(equations.gamma * p_rr / rho_rr)
    else # y-direction
        λ_min = v2_ll - sqrt(equations.gamma * p_ll / rho_ll)
        λ_max = v2_rr + sqrt(equations.gamma * p_rr / rho_rr)
    end

    return λ_min, λ_max
end

@inline function min_max_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    norm_ = norm(normal_direction)
    # The v_normals are already scaled by the norm
    λ_min = v_normal_ll - sqrt(equations.gamma * p_ll / rho_ll) * norm_
    λ_max = v_normal_rr + sqrt(equations.gamma * p_rr / rho_rr) * norm_

    return λ_min, λ_max
end

@inline function max_abs_speeds(u,
                                equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, v1, v2, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c
end

# Convert conservative variables to primitive
@inline function cons2prim(u,
                           equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e, phi = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2))

    return SVector(rho, v1, v2, p, phi)
end

# Convert conservative variables to entropy
@inline function cons2entropy(u,
                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v_square = v1^2 + v2^2
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * v_square)
    s = log(p) - equations.gamma * log(rho)
    rho_p = rho / p

    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         0.5f0 * rho_p * v_square
    w2 = rho_p * v1
    w3 = rho_p * v2
    w4 = -rho_p

    return SVector(w1, w2, w3, w4, zero(eltype(u)))
end

@inline function entropy2cons(w,
                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # See Hughes, Franca, Mallet (1986) A new finite element formulation for CFD
    # [DOI: 10.1016/0045-7825(86)90127-1](https://doi.org/10.1016/0045-7825(86)90127-1)
    @unpack gamma = equations

    # convert to entropy `-rho * s` used by Hughes, France, Mallet (1986)
    # instead of `-rho * s / (gamma - 1)`
    V1, V2, V3, V5 = w .* (gamma - 1)

    # s = specific entropy, eq. (53)
    s = gamma - V1 + (V2^2 + V3^2) / (2 * V5)

    # eq. (52)
    rho_iota = ((gamma - 1) / (-V5)^gamma)^(equations.inv_gamma_minus_one) *
               exp(-s * equations.inv_gamma_minus_one)

    # eq. (51)
    rho = -rho_iota * V5
    rho_v1 = rho_iota * V2
    rho_v2 = rho_iota * V3
    rho_e = rho_iota * (1 - (V2^2 + V3^2) / (2 * V5))
    return SVector(rho, rho_v1, rho_v2, rho_e, zero(eltype(w)))
end

# Convert primitive to conservative variables
@inline function prim2cons(prim,
                           equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, v1, v2, p, phi = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_e = p * equations.inv_gamma_minus_one + 0.5f0 * (rho_v1 * v1 + rho_v2 * v2)
    return SVector(rho, rho_v1, rho_v2, rho_e, phi)
end

@inline function pressure(u,
                          equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e = u
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1^2 + rho_v2^2) / rho)
    return p
end

# Calculate thermodynamic entropy for a conservative state `cons`
@inline function entropy_thermodynamic(cons,
                                       equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Pressure
    p = (equations.gamma - 1) * (cons[4] - 0.5f0 * (cons[2]^2 + cons[3]^2) / cons[1])

    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])

    return s
end

# Calculate mathematical entropy for a conservative state `cons`
@inline function entropy_math(cons,
                              equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    # Mathematical entropy
    S = -entropy_thermodynamic(cons, equations) * cons[1] *
        equations.inv_gamma_minus_one

    return S
end

# Default entropy is the mathematical entropy
@inline function entropy(cons,
                         equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    entropy_math(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline energy_total(cons, ::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D) = cons[4]

# Calculate kinetic energy for a conservative state `cons`
@inline function energy_kinetic(u,
                                equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e = u
    return (rho_v1^2 + rho_v2^2) / (2 * rho)
end

# Calculate internal energy for a conservative state `cons`
@inline function energy_internal(cons,
                                 equations::CompressibleEulerInternalKineticEnergyEquationsWithGravity2D)
    return energy_total(cons, equations) - energy_kinetic(cons, equations)
end
end # @muladd
