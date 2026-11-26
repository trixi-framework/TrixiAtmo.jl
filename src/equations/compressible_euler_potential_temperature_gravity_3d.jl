@muladd begin
#! format: noindent

@doc raw"""
    CompressibleEulerPotentialTemperatureEquationsWithGravity3D(; c_p, c_v, gravity,
                                                                turbulence = NoTurbulence())

The compressible Euler equations with gravity in potential temperature formulation
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho v_2 \\ \rho v_3 \\ \rho \theta
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
\rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ \rho v_1 v_3 \\ \rho \theta v_1
\end{pmatrix}
+
\frac{\partial}{\partial y}
\begin{pmatrix}
\rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ \rho v_2 v_3 \\ \rho \theta v_2
\end{pmatrix}
+
\frac{\partial}{\partial z}
\begin{pmatrix}
\rho v_3 \\ \rho v_1 v_3 \\ \rho v_2 v_3 \\ \rho v_3^2 + p \\ \rho \theta v_3
\end{pmatrix}
=
\begin{pmatrix}
0 \\ -\rho \frac{\partial}{\partial x} \phi \\ -\rho \frac{\partial}{\partial y} \phi \\ -\rho \frac{\partial}{\partial z} \phi \\ 0
\end{pmatrix}
```
for an ideal gas with ratio of specific heats ``\gamma = c_p / c_v`` in three space dimensions.
Here, ``\rho`` is the density, ``v_1, v_2, v_3`` are the velocities, ``\theta`` is the
potential temperature, ``\phi`` is the gravitational potential, and
```math
p = p_0 \left( \frac{R \rho \theta}{p_0} \right)^\gamma
```
the pressure, where ``p_0 = 10^5\,\mathrm{Pa}`` is the reference pressure and
``R = c_p - c_v`` the gas constant.

The geopotential ``\phi`` is the only auxiliary variable, see `Trixi.nauxvariables`. Its
values are given by the function `aux_field(x, equations)` passed to
`SemidiscretizationHyperbolic`, e.g., [`aux_field_gravity_earth`](@ref). Thus,
the pointwise functions take the auxiliary variables `aux = SVector(phi)` as additional
argument after the state, e.g., `flux(u, aux, normal_direction, equations)` or
`flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)`.

The turbulent fluxes of the model `turbulence`, e.g., a
[`HorizontalVerticalTurbulence`](@ref), are the parabolic terms of
[`TurbulentDiffusion3D`](@ref). The prognostic variables of the turbulence model, e.g., the
density of the turbulent kinetic energy of the [`DeardorffTKE`](@ref), follow
``\rho \theta`` and are advected as specific quantities.
"""
struct CompressibleEulerPotentialTemperatureEquationsWithGravity3D{NVARS, RealT <: Real,
                                                                   Turbulence} <:
       AbstractCompressibleEulerEquations{3, NVARS}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
    turbulence::Turbulence
    function CompressibleEulerPotentialTemperatureEquationsWithGravity3D(; c_p, c_v,
                                                                         gravity,
                                                                         turbulence = NoTurbulence())
        c_p, c_v, g = promote(c_p, c_v, gravity)
        turbulence = similar(turbulence, typeof(c_p))
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        K = p_0 * (R / p_0)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        NVARS = 5 + nturbulencevariables(turbulence)
        return new{NVARS, typeof(c_p), typeof(turbulence)}(p_0, c_p, c_v, g, R,
                                                           gamma,
                                                           inv_gamma_minus_one,
                                                           K, stolarsky_factor,
                                                           turbulence)
    end
end

# Together with the specialization of `Adapt.adapt_structure` in Trixi.jl,
# this allows to move semidiscretizations and their components including
# the equations to GPUs and adapt the floating point type, e.g.,
# to `Float32` to improve performance on GPUs.
function Base.similar(equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
                      ::Type{NewRealT}) where {NewRealT}
    c_p, c_v, gravity = convert.(NewRealT, (equations.c_p, equations.c_v, equations.g))
    turbulence = similar(equations.turbulence, NewRealT)
    return CompressibleEulerPotentialTemperatureEquationsWithGravity3D(; c_p, c_v,
                                                                       gravity,
                                                                       turbulence)
end

function varnames(::typeof(cons2cons),
                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    ("rho", "rho_v1", "rho_v2", "rho_v3", "rho_theta",
     varnames_turbulence(cons2cons, equations.turbulence)...)
end

function varnames(::typeof(cons2prim),
                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    ("rho", "v1", "v2", "v3", "p1",
     varnames_turbulence(cons2prim, equations.turbulence)...)
end

@inline Trixi.nauxvariables(::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = 1

Trixi.auxvarnames(::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = ("phi",)

have_nonconservative_terms(::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = Trixi.True()

# The variables of the turbulence model follow the density, the momentum, and `rho_theta`
@inline function turbulence_variables(u,
                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    nturbulence = nvariables(equations) - 5
    return SVector{nturbulence, eltype(u)}(ntuple(i -> u[5 + i], Val(nturbulence)))
end

# Flux of the variables of the turbulence model advected as specific quantities with the
# arithmetic mean, for the mass flux `f1` of a two-point flux
@inline function flux_turbulence_central(f1, u_ll, u_rr,
                                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    return f1 * 0.5f0 *
           (turbulence_variables(u_ll, equations) / u_ll[1] +
            turbulence_variables(u_rr, equations) / u_rr[1])
end

# There are no separate auxiliary variables outside of the domain, so the inner ones are
# used on both sides.
@inline function boundary_condition_slip_wall(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_functions,
                                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4] - 2 * u_normal * normal[3],
                         u_inner[5], turbulence_variables(u_inner, equations)...)

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, aux_inner, aux_inner,
                                 normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary,
                                                 aux_inner, aux_inner,
                                                 normal_direction, equations)
    return flux, noncons_flux
end

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, aux, normal_direction::AbstractVector,
                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
               v3 * normal_direction[3]
    rho_v_normal = rho * v_normal
    p = pressure(u, aux, equations)

    f1 = rho_v_normal
    f2 = rho_v_normal * v1 + p * normal_direction[1]
    f3 = rho_v_normal * v2 + p * normal_direction[2]
    f4 = rho_v_normal * v3 + p * normal_direction[3]
    f5 = rho_theta * v_normal
    f_turbulence = turbulence_variables(u, equations) * v_normal
    return SVector(f1, f2, f3, f4, f5, f_turbulence...)
end

"""
    flux_nonconservative_isothermal(u_ll, u_rr, aux_ll, aux_rr,
                                    normal_direction::AbstractVector,
                                    equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Well-balanced gravity term for an isothermal background state
for the [`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref)
developed by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_isothermal(u_ll, u_rr, aux_ll, aux_rr,
                                                 normal_direction::AbstractVector,
                                                 equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll = u_ll[1]
    rho_rr = u_rr[1]
    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    rho_avg = ln_mean(rho_ll, rho_rr)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   normal_direction[3] * rho_avg * phi_jump,
                   zero(eltype(u_ll)),
                   zero(turbulence_variables(u_ll, equations))...)
end

"""
    flux_nonconservative_artiano_etal(u_ll, u_rr, aux_ll, aux_rr,
                                      normal_direction::AbstractVector,
                                      equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Well-balanced gravity term for constant potential temperature background state by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_artiano_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                   normal_direction::AbstractVector,
                                                   equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll = u_ll[1]
    rho_rr = u_rr[1]
    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   normal_direction[3] * rho_avg * phi_jump,
                   zero(eltype(u_ll)),
                   zero(turbulence_variables(u_ll, equations))...)
end

"""
    flux_nonconservative_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                    normal_direction::AbstractVector,
                                    equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

-  Souza et al.
   The Flux-Differencing Discontinuous {G}alerkin Method Applied to
   an Idealized Fully Compressible Nonhydrostatic Dry Atmosphere
   [DOI: 10.1029/2022MS003527] (https://doi.org/10.1029/2022MS003527)
"""
@inline function flux_nonconservative_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                 normal_direction::AbstractVector,
                                                 equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll = u_ll[1]
    rho_rr = u_rr[1]
    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    phi_jump = phi_rr - phi_ll
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   normal_direction[3] * rho_avg * phi_jump,
                   zero(eltype(u_ll)),
                   zero(turbulence_variables(u_ll, equations))...)
end

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, aux_ll, aux_rr,
                                         normal_direction::AbstractVector,
                                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    a = flux_lmars.speed_of_sound
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f = u_ll * v_interface
    else
        f = u_rr * v_interface
    end

    return SVector(f[1],
                   f[2] + p_interface * normal_direction[1],
                   f[3] + p_interface * normal_direction[2],
                   f[4] + p_interface * normal_direction[3],
                   f[5], turbulence_variables(f, equations)...)
end

"""
    flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                        equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Version of the low Mach number approximate Riemann solver `FluxLMARS(340)` that returns the
conservative and the nonconservative contribution of the surface flux as a tuple, see
[`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref). It is therefore used as a
single `surface_flux` instead of the tuple `(FluxLMARS(340), flux_zero)`.

The geopotential is continuous across interfaces, so the nonconservative gravity term, which
is proportional to the jump of the geopotential, vanishes on the surface. Hence both
contributions are given by the plain LMARS flux and no gravity term is computed here.
"""
@inline function flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # constants
    a = 340

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    norm_ = norm(normal_direction)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho_avg * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho_avg) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f = u_ll * v_interface
    else
        f = u_rr * v_interface
    end

    flux = SVector(f[1],
                   f[2] + p_interface * normal_direction[1],
                   f[3] + p_interface * normal_direction[2],
                   f[4] + p_interface * normal_direction[3],
                   f[5], turbulence_variables(f, equations)...)

    return flux, flux
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_lmars_combined), ::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = Trixi.True()

"""
    flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Entropy-stable and well-balanced version of the low Mach number approximate Riemann solver
[`FluxLMARS`](@ref) by
- Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
  Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
  Euler Equations in Vector-Invariant Form
  [arXiv: 2610.05820](https://arxiv.org/abs/2610.05820)

Like [`flux_lmars_combined`](@ref), it returns the conservative and the nonconservative
contribution of the surface flux as a tuple, see
[`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref).
"""
@inline function flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr,
                               normal_direction::AbstractVector,
                               equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)
    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]

    norm_ = norm(normal_direction)

    # Arithmetic mean of the local speeds of sound instead of a constant
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)
    c_avg = 0.5f0 * (c_ll + c_rr)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    rho_mean = ln_mean(rho_ll, rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) -
                  0.5f0 * c_avg * rho_avg * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) -
                  0.5f0 / (c_avg * rho_mean) * (p_rr - p_ll) * norm_

    # Upwinded mass, momentum, and turbulence fluxes
    if v_interface > 0
        f1 = u_ll[1] * v_interface
        f2 = u_ll[2] * v_interface
        f3 = u_ll[3] * v_interface
        f4 = u_ll[4] * v_interface
        f_turbulence = turbulence_variables(u_ll, equations) * v_interface
        sign_v = one(f1)
    else
        f1 = u_rr[1] * v_interface
        f2 = u_rr[2] * v_interface
        f3 = u_rr[3] * v_interface
        f4 = u_rr[4] * v_interface
        f_turbulence = turbulence_variables(u_rr, equations) * v_interface
        sign_v = -one(f1)
    end

    # Entropy-stable flux of the potential temperature: the entropy-conservative mean
    # `1 / {{1/theta}}_ln` of `flux_ec` and an upwind correction scaled by the jump
    theta_ll = rho_theta_ll / rho_ll
    theta_rr = rho_theta_rr / rho_rr
    theta_mean = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    theta_min = min(theta_ll, theta_rr)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    f5 = f1 *
         (theta_mean - theta_min / (2 * theta_avg) * (theta_rr - theta_ll) * sign_v)

    flux = SVector(f1,
                   f2 + p_interface * normal_direction[1],
                   f3 + p_interface * normal_direction[2],
                   f4 + p_interface * normal_direction[3],
                   f5, f_turbulence...)

    return flux, flux
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_lmars_es), ::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = Trixi.True()

# The boundary condition must also be specialized for a combined surface flux such that it
# returns only the first contribution of the surface flux.
@inline function boundary_condition_slip_wall(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function::Union{typeof(flux_lmars_combined),
                                                                           typeof(flux_lmars_es)},
                                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4] - 2 * u_normal * normal[3],
                         u_inner[5], turbulence_variables(u_inner, equations)...)

    # calculate the boundary flux
    flux, _ = surface_flux_function(u_inner, u_boundary, aux_inner, aux_inner,
                                    normal_direction, equations)

    return flux
end

"""
    flux_kennedy_gruber(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                        equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Modification of the kinetic-energy-preserving two-point flux by Kennedy and Gruber
that is consistent with [`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref). See:
- Kennedy and Gruber (2008)
  Reduced aliasing formulations of the convective terms within the
  Navier-Stokes equations for a compressible fluid
  [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function flux_kennedy_gruber(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Unpack left and right state
    rho_theta_ll = u_ll[5]
    rho_theta_rr = u_rr[5]
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    theta_avg = 0.5f0 * (rho_theta_ll / rho_ll + rho_theta_rr / rho_rr)

    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2] +
                  v3_avg * normal_direction[3]

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * theta_avg
    f_turbulence = flux_turbulence_central(f1, u_ll, u_rr, equations)

    return SVector(f1, f2, f3, f4, f5, f_turbulence...)
end

"""
    flux_kennedy_gruber_souza_etal(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                                   equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Fused version of [`flux_kennedy_gruber`](@ref) and [`flux_nonconservative_souza_etal`](@ref)
that computes the conservative and the nonconservative contributions in a single kernel,
returning the tuple

    flux_kennedy_gruber(u_ll, u_rr, aux_ll, aux_rr, n, equations)
        ± 0.5f0 * flux_nonconservative_souza_etal(u_ll, u_rr, aux_ll, aux_rr, n, equations)

see [`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref). It is therefore used as
a single `volume_flux` instead of the tuple
`(flux_kennedy_gruber, flux_nonconservative_souza_etal)`.
"""
@inline function flux_kennedy_gruber_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                normal_direction::AbstractVector,
                                                equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Unpack left and right state
    rho_theta_ll = u_ll[5]
    rho_theta_rr = u_rr[5]
    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    theta_avg = 0.5f0 * (rho_theta_ll / rho_ll + rho_theta_rr / rho_rr)

    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2] +
                  v3_avg * normal_direction[3]

    # Conservative part, cf. `flux_kennedy_gruber`
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * theta_avg

    # Nonconservative part, cf. `flux_nonconservative_souza_etal`, scaled by 0.5
    gravity = 0.5f0 * rho_avg * (phi_rr - phi_ll)
    g2 = normal_direction[1] * gravity
    g3 = normal_direction[2] * gravity
    g4 = normal_direction[3] * gravity

    f_turbulence = flux_turbulence_central(f1, u_ll, u_rr, equations)

    return SVector(f1, f2 + g2, f3 + g3, f4 + g4, f5, f_turbulence...),
           SVector(f1, f2 - g2, f3 - g3, f4 - g4, f5, f_turbulence...)
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_kennedy_gruber_souza_etal), ::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = Trixi.True()

"""
    flux_tec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
             equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_tec(u_ll, u_rr, aux_ll, aux_rr,
                          normal_direction::AbstractVector,
                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]

    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f_turbulence = flux_turbulence_central(f1, u_ll, u_rr, equations)
    return SVector(f1, f2, f3, f4, f5, f_turbulence...)
end

"""
    flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
            equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Entropy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_ec(u_ll, u_rr, aux_ll, aux_rr,
                         normal_direction::AbstractVector,
                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]
    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr
    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    f_turbulence = flux_turbulence_central(f1, u_ll, u_rr, equations)
    return SVector(f1, f2, f3, f4, f5, f_turbulence...)
end

"""
    flux_ec_isothermal(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
                       equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Fused version of [`flux_ec`](@ref) and [`flux_nonconservative_isothermal`](@ref) that
returns the conservative and the nonconservative contributions as the tuple

    flux_ec(u_ll, u_rr, aux_ll, aux_rr, n, equations)
        ± 0.5f0 * flux_nonconservative_isothermal(u_ll, u_rr, aux_ll, aux_rr, n, equations)

see [`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref). It is therefore used as
a single `volume_flux` instead of the tuple `(flux_ec, flux_nonconservative_isothermal)`,
e.g., on GPUs.
"""
@inline function flux_ec_isothermal(u_ll, u_rr, aux_ll, aux_rr,
                                    normal_direction::AbstractVector,
                                    equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    flux = flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)
    # Nonconservative part, cf. `flux_nonconservative_isothermal`, scaled by 0.5
    gravity = 0.5f0 * ln_mean(u_ll[1], u_rr[1]) * (aux_rr[1] - aux_ll[1])
    zero_ = zero(gravity)
    g = SVector(zero_, normal_direction[1] * gravity, normal_direction[2] * gravity,
                normal_direction[3] * gravity, zero_,
                zero(turbulence_variables(u_ll, equations))...)
    return flux + g, flux - g
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_ec_isothermal), ::CompressibleEulerPotentialTemperatureEquationsWithGravity3D) = Trixi.True()

"""
    flux_etec(u_ll, u_rr, aux_ll, aux_rr, normal_direction,
              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)

Entropy and total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_etec(u_ll, u_rr, aux_ll, aux_rr,
                           normal_direction::AbstractVector,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]
    _, _, _, _, rho_theta_ll = u_ll
    _, _, _, _, rho_theta_rr = u_rr
    # Compute the necessary mean values

    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f1 = f5 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f_turbulence = flux_turbulence_central(f1, u_ll, u_rr, equations)

    return SVector(f1, f2, f3, f4, f5, f_turbulence...)
end

@inline function Trixi.prim2cons(prim, aux,
                                 equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho, v1, v2, v3, p = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_v3 = rho * v3
    rho_theta = (p / equations.K)^(1 / equations.gamma)
    return SVector(rho, rho_v1, rho_v2, rho_v3, rho_theta,
                   (rho * turbulence_variables(prim, equations))...)
end

@inline function cons2prim(u, aux,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = equations.K * rho_theta^equations.gamma
    return SVector(rho, v1, v2, v3, p, (turbulence_variables(u, equations) / rho)...)
end

@inline function cons2cons(u, aux,
                           equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    return u
end

@inline function cons2entropy(u, aux,
                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta = u

    w1 = log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma
    w5 = rho / rho_theta * equations.gamma

    return SVector(w1, zero(eltype(u)), zero(eltype(u)), zero(eltype(u)), w5,
                   zero(turbulence_variables(u, equations))...)
end

@inline function entropy(cons, aux,
                         equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    p = equations.K * cons[5]^equations.gamma
    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])
    S = -s * cons[1] / (equations.gamma - 1)
    return S
end

@inline function pressure(cons, aux,
                          equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    p = equations.K * exp(equations.gamma * log(cons[5]))
    return p
end

# Calculate kinetic energy for a conservative state `cons`
@inline function energy_kinetic(u, aux,
                                equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho, rho_v1, rho_v2, rho_v3, _ = u
    return 0.5f0 * (rho_v1^2 + rho_v2^2 + rho_v3^2) / rho
end

@inline function energy_total(cons, aux,
                              equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Mathematical entropy
    p = equations.p_0 * (equations.R * cons[5] / equations.p_0)^equations.gamma

    U = (p / (equations.gamma - 1) +
         0.5f0 * (cons[2]^2 + cons[3]^2 + cons[4]^2) / (cons[1]))

    return U
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr, orientation::Integer,
                                     equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    else # orientation == 3
        v_ll = v3_ll
        v_rr = v3_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1]
            + v2_ll * normal_direction[2]
            + v3_ll * normal_direction[3])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            + v2_rr * normal_direction[2]
            + v3_rr * normal_direction[3])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, aux_ll, aux_rr, orientation::Integer,
                               equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    else # orientation == 3
        v_ll = v3_ll
        v_rr = v3_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, aux_ll, aux_rr,
                               normal_direction::AbstractVector,
                               equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Calculate normal velocities and sound speeds
    # left
    v_ll = (v1_ll * normal_direction[1]
            + v2_ll * normal_direction[2]
            + v3_ll * normal_direction[3])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
            + v2_rr * normal_direction[2]
            + v3_rr * normal_direction[3])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_, abs(v_rr) + c_rr * norm_)
end

@inline function max_abs_speeds(u, aux,
                                equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    rho, v1, v2, v3, p = cons2prim(u, aux, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c, abs(v3) + c
end
end # @muladd

# Specializations of the `FluxTurbo` interface of Trixi.jl for the volume fluxes
# `FluxTurbo(flux_kennedy_gruber, flux_nonconservative_souza_etal)` and
# `FluxTurbo(flux_ec, flux_nonconservative_isothermal)` without variables of a turbulence
# model, see below for the density of the TKE.
# We do not wrap this code in `@muladd begin ... end` block. Optimizations like
# this are handled automatically by LoopVectorization.jl.
const FluxKennedyGruberSouzaEtal = Tuple{typeof(flux_kennedy_gruber),
                                         typeof(flux_nonconservative_souza_etal)}
const FluxECIsothermal = Tuple{typeof(flux_ec),
                               typeof(flux_nonconservative_isothermal)}

# We precompute `rho, v1, v2, v3, theta, p, phi`, i.e., one quantity more than the number
# of variables and auxiliary variables. Both fluxes use the same precomputed variables.
@inline function Trixi.nturbovars(::Union{FluxKennedyGruberSouzaEtal,
                                          FluxECIsothermal},
                                  ::CompressibleEulerPotentialTemperatureEquationsWithGravity3D{5})
    return Val(7)
end

@inline function Trixi.cons2turbo(::Union{FluxKennedyGruberSouzaEtal,
                                          FluxECIsothermal},
                                  rho, rho_v1, rho_v2, rho_v3, rho_theta, phi,
                                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D{5})
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    theta = rho_theta / rho
    p = equations.K * exp(equations.gamma * log(rho_theta))

    return rho, v1, v2, v3, theta, p, phi
end

@inline function flux_turbo_dry_air(::FluxKennedyGruberSouzaEtal,
                                    rho_ll, v1_ll, v2_ll, v3_ll, theta_ll, p_ll, phi_ll,
                                    rho_rr, v1_rr, v2_rr, v3_rr, theta_rr, p_rr, phi_rr,
                                    normal_direction_1, normal_direction_2,
                                    normal_direction_3,
                                    equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)

    v_dot_n_avg = v1_avg * normal_direction_1 + v2_avg * normal_direction_2 +
                  v3_avg * normal_direction_3

    # Conservative part, i.e. 'flux_kennedy_gruber'
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction_1
    f3 = f1 * v2_avg + p_avg * normal_direction_2
    f4 = f1 * v3_avg + p_avg * normal_direction_3
    f5 = f1 * theta_avg

    # Nonconservative part, i.e., `flux_nonconservative_souza_etal`, scaled by 0.5. Since this
    # flux is antisymmetric, the right contribution is obtained by flipping the sign.
    gravity = 0.5f0 * rho_avg * (phi_rr - phi_ll)
    g2 = normal_direction_1 * gravity
    g3 = normal_direction_2 * gravity
    g4 = normal_direction_3 * gravity

    flux_left = SVector(f1, f2 + g2, f3 + g3, f4 + g4, f5)
    flux_right = SVector(f1, f2 - g2, f3 - g3, f4 - g4, f5)

    return flux_left, flux_right
end

@inline function Trixi.flux_turbo(volume_flux::Union{FluxKennedyGruberSouzaEtal,
                                                     FluxECIsothermal},
                                  rho_ll, v1_ll, v2_ll, v3_ll, theta_ll, p_ll, phi_ll,
                                  rho_rr, v1_rr, v2_rr, v3_rr, theta_rr, p_rr, phi_rr,
                                  normal_direction_1, normal_direction_2,
                                  normal_direction_3,
                                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D{5})
    return flux_turbo_dry_air(volume_flux,
                              rho_ll, v1_ll, v2_ll, v3_ll, theta_ll, p_ll, phi_ll,
                              rho_rr, v1_rr, v2_rr, v3_rr, theta_rr, p_rr, phi_rr,
                              normal_direction_1, normal_direction_2, normal_direction_3,
                              equations)
end

# Fused version of `flux_ec` and `flux_nonconservative_isothermal`, see
# - Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
#   Structure-Preserving High-Order Methods for the Compressible Euler Equations
#   in Potential Temperature Formulation for Atmospheric Flows
#   (https://arxiv.org/abs/2509.10311)
@inline function flux_turbo_dry_air(::FluxECIsothermal,
                                    rho_ll, v1_ll, v2_ll, v3_ll, theta_ll, p_ll, phi_ll,
                                    rho_rr, v1_rr, v2_rr, v3_rr, theta_rr, p_rr, phi_rr,
                                    normal_direction_1, normal_direction_2,
                                    normal_direction_3,
                                    equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
    v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
                 v3_ll * normal_direction_3
    v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
                 v3_rr * normal_direction_3

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Conservative part, i.e., `flux_ec`. Since `rho / rho_theta = 1 / theta`, the
    # potential temperature flux uses the same mean value as in `flux_ec`.
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction_1
    f3 = f1 * v2_avg + p_avg * normal_direction_2
    f4 = f1 * v3_avg + p_avg * normal_direction_3
    f5 = f1 * inv_ln_mean(1 / theta_ll, 1 / theta_rr)

    # Since this flux is antisymmetric, the right contribution is obtained by flipping the
    # sign.
    gravity = 0.5f0 * rho_mean * (phi_rr - phi_ll)
    g2 = normal_direction_1 * gravity
    g3 = normal_direction_2 * gravity
    g4 = normal_direction_3 * gravity

    flux_left = SVector(f1, f2 + g2, f3 + g3, f4 + g4, f5)
    flux_right = SVector(f1, f2 - g2, f3 - g3, f4 - g4, f5)

    return flux_left, flux_right
end

# Specializations of the `FluxTurbo` interface for the equations with the density of the
# turbulent kinetic energy `rho_e` of the `DeardorffTKE` as sixth variable: the precomputed
# variables are the ones without it followed by `e = rho_e / rho`, which is advected with the
# mass flux and the arithmetic mean as in `flux_turbulence_central`
@inline function Trixi.nturbovars(::Union{FluxKennedyGruberSouzaEtal,
                                          FluxECIsothermal},
                                  ::CompressibleEulerPotentialTemperatureEquationsWithGravity3D{6})
    return Val(8)
end

@inline function Trixi.cons2turbo(volume_flux::Union{FluxKennedyGruberSouzaEtal,
                                                     FluxECIsothermal},
                                  rho, rho_v1, rho_v2, rho_v3, rho_theta, rho_e, phi,
                                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D{6})
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    theta = rho_theta / rho
    p = equations.K * exp(equations.gamma * log(rho_theta))

    return rho, v1, v2, v3, theta, p, phi, rho_e / rho
end

@inline function Trixi.flux_turbo(volume_flux::Union{FluxKennedyGruberSouzaEtal,
                                                     FluxECIsothermal},
                                  rho_ll, v1_ll, v2_ll, v3_ll, theta_ll, p_ll, phi_ll, e_ll,
                                  rho_rr, v1_rr, v2_rr, v3_rr, theta_rr, p_rr, phi_rr, e_rr,
                                  normal_direction_1, normal_direction_2,
                                  normal_direction_3,
                                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D{6})
    flux_left, flux_right = flux_turbo_dry_air(volume_flux,
                                               rho_ll, v1_ll, v2_ll, v3_ll, theta_ll,
                                               p_ll, phi_ll,
                                               rho_rr, v1_rr, v2_rr, v3_rr, theta_rr,
                                               p_rr, phi_rr,
                                               normal_direction_1, normal_direction_2,
                                               normal_direction_3, equations)
    f_e = flux_left[1] * 0.5f0 * (e_ll + e_rr)
    return SVector(flux_left..., f_e), SVector(flux_right..., f_e)
end
