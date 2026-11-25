@muladd begin
#! format: noindent

@doc raw"""
    CompressibleEulerMoistPotentialTemperatureEquations3D(; c_p, c_v, gravity,
                                                          R_v = 461.5,
                                                          microphysics = KesslerMicrophysics())

The compressible Euler equations with gravity in potential temperature formulation with
moisture, whose moisture variables, precipitation, and right-hand sides are given by the
microphysics scheme `microphysics`, e.g., [`KesslerMicrophysics`](@ref) with water vapor,
cloud water, and rain water as in the test cases of DCMIP2016,
```math
\begin{aligned}
\partial_t \varrho + \nabla \cdot (\varrho \vec{v} - F \hat{k}) &= 0, \\
\partial_t (\varrho \vec{v}) + \nabla \cdot (\varrho \vec{v} \otimes \vec{v} + p I
  - F \vec{v} \otimes \hat{k}) &= - \varrho \nabla \phi, \\
\partial_t (\varrho \theta_v) + \nabla \cdot (\varrho \theta_v \vec{v})
  &= \varrho_d \left(1 + \frac{R_v}{R} q_v\right) \frac{\partial_t T}{\pi}
  + \frac{R_v}{R} \theta S_v, \\
\partial_t \varrho_x + \nabla \cdot (\varrho_x \vec{v} - F_x \hat{k}) &= S_x,
\end{aligned}
```
for each moisture variable ``\varrho_x``, the first of which is the density of water vapor
``\varrho_v``. Here, the total density ``\varrho = \varrho_d + \sum_x \varrho_x`` sums the
density of dry air and the moisture variables that are densities of water, ``F_x`` is the
downward flux of the precipitation of ``\varrho_x`` through a horizontal surface, ``F`` the
sum of the ones of the densities of water, and ``S_x`` and the isobaric heating rate
``\partial_t T`` are the right-hand sides of the microphysics, see
[`source_terms_microphysics`](@ref). The virtual potential temperature
```math
\theta_v = \frac{p}{\varrho R \pi} = \theta \frac{1 + q_v R_v / R}{1 + \sum_x q_x}
```
includes the weight of the condensed water, with the gas constants ``R = c_p - c_v`` of dry
air and ``R_v`` of water vapor, the dry potential temperature ``\theta = T / \pi``, the
Exner pressure ``\pi = (p / p_0)^{R / c_p}``, and the mixing ratios
``q_x = \varrho_x / \varrho_d``. The precipitation falls in the direction ``-\hat{k}``,
where ``\hat{k}`` is the upward unit vector. The pressure follows from Dalton's law
``p = (\varrho_d R + \varrho_v R_v) T``, i.e., ``p = K (\varrho \theta_v)^\gamma``, which
depends on the fifth variable (``\theta_v``) only.

The auxiliary variables, see `Trixi.nauxvariables`, are the geopotential ``\phi`` and the
Cartesian components of ``\hat{k}``, i.e., `aux = SVector(phi, k1, k2, k3)`. Their values
are given by the function `aux_field(x, equations)` passed to
`SemidiscretizationHyperbolic`. On the sphere, ``\hat{k} = \vec{x} / \lVert \vec{x} \rVert``
is the radial direction, which varies within each element. The pointwise functions take the
auxiliary variables as additional argument after the state, e.g.,
`flux(u, aux, normal_direction, equations)` and `cons2prim(u, aux, equations)`.
The numerical fluxes return the conservative and the nonconservative contribution together,
see [`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref), i.e., use
[`flux_kennedy_gruber_souza_etal`](@ref) as `volume_flux` and
[`flux_lmars_combined`](@ref) as `surface_flux`.
"""
struct CompressibleEulerMoistPotentialTemperatureEquations3D{NVARS, RealT <: Real,
                                                             Microphysics} <:
       AbstractCompressibleEulerEquations{3, NVARS}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat of dry air at constant pressure in J/(kg K)
    c_v::RealT # specific heat of dry air at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant of dry air in J/(kg K)
    R_v::RealT # gas constant of water vapor in J/(kg K)
    gamma::RealT # ratio of specific heats
    # Scaling factor between the pressure and the weighted potential temperature
    K::RealT # = p_0 * (R / p_0)^gamma
    microphysics::Microphysics
    function CompressibleEulerMoistPotentialTemperatureEquations3D(; c_p, c_v, gravity,
                                                                   R_v = 461.5,
                                                                   microphysics = KesslerMicrophysics())
        c_p, c_v, g = promote(c_p, c_v, gravity)
        RealT = typeof(c_p)
        R_v = convert(RealT, R_v)
        microphysics = similar(microphysics, RealT)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        K = p_0 * (R / p_0)^gamma
        NVARS = 5 + nmoisturevariables(microphysics)
        return new{NVARS, RealT, typeof(microphysics)}(p_0, c_p, c_v, g, R, R_v, gamma,
                                                       K,
                                                       microphysics)
    end
end

# Together with the specialization of `Adapt.adapt_structure` in Trixi.jl,
# this allows to move semidiscretizations and their components including
# the equations to GPUs and adapt the floating point type, e.g.,
# to `Float32` to improve performance on GPUs.
function Base.similar(equations::CompressibleEulerMoistPotentialTemperatureEquations3D,
                      ::Type{NewRealT}) where {NewRealT}
    c_p, c_v, gravity, R_v = convert.(NewRealT,
                                      (equations.c_p, equations.c_v, equations.g,
                                       equations.R_v))
    microphysics = similar(equations.microphysics, NewRealT)
    return CompressibleEulerMoistPotentialTemperatureEquations3D(; c_p, c_v, gravity,
                                                                 R_v,
                                                                 microphysics)
end

function varnames(::typeof(cons2cons),
                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return ("rho", "rho_v1", "rho_v2", "rho_v3", "rho_theta_v",
            varnames_moisture(cons2cons, equations.microphysics)...)
end

# The primitive variables use the mixing ratios with respect to dry air, as the Kessler
# microphysics of DCMIP2016
function varnames(::typeof(cons2prim),
                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return ("rho", "v1", "v2", "v3", "p",
            varnames_moisture(cons2prim, equations.microphysics)...)
end

@inline Trixi.nauxvariables(::CompressibleEulerMoistPotentialTemperatureEquations3D) = 4

Trixi.auxvarnames(::CompressibleEulerMoistPotentialTemperatureEquations3D) = ("phi",
                                                                              "k1",
                                                                              "k2",
                                                                              "k3")

have_nonconservative_terms(::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# The moisture variables follow the density, the momentum, and `rho_theta_v`
@inline function moisture_variables(u,
                                    equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return SVector(ntuple(i -> u[5 + i], Val(nvariables(equations) - 5)))
end

# Slip-wall boundary condition
# As usual, the external state mirrors the normal velocity, such that there is no flux of
# mass through the boundary. In addition, the precipitation leaves the domain through the
# lower boundary with the sedimentation flux of the inner state, while no precipitation
# falls into the domain through the upper boundary. There are no separate auxiliary
# variables outside of the domain, so the inner ones are used on both sides.
@inline function boundary_condition_slip_wall(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function::Union{typeof(flux_lmars_combined),
                                                                           typeof(flux_lmars_es)},
                                              equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal momentum
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4] - 2 * u_normal * normal[3],
                         u_inner[5], moisture_variables(u_inner, equations)...)

    # calculate the boundary flux
    flux_without_sedimentation = surface_flux_without_sedimentation(surface_flux_function)
    flux = flux_without_sedimentation(u_inner, u_boundary, aux_inner, aux_inner,
                                      normal_direction, equations)
    if vertical_component(aux_inner, normal_direction, equations) < 0
        flux = flux -
               flux_sedimentation(u_inner, aux_inner, normal_direction, equations)
    end

    return flux
end

# Component of the `normal_direction` in the direction of the upward unit vector `k`, which
# is given by the auxiliary variables `aux[2:4]` in the physical coordinates
@inline function vertical_component(aux, normal_direction::AbstractVector,
                                    equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return aux[2] * normal_direction[1] + aux[3] * normal_direction[2] +
           aux[4] * normal_direction[3]
end

# Flux of the precipitation through a face with the normal direction `normal_direction`,
# i.e., of the total density, of the momentum carried by the precipitation with the
# velocity of the state `u`, and of the moisture variables
@inline function flux_sedimentation(u, aux, normal_direction::AbstractVector,
                                    equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho = u[1]
    precipitation, moisture = sedimentation(moisture_variables(u, equations),
                                            density_dry(u, equations),
                                            equations.microphysics)
    k_normal = vertical_component(aux, normal_direction, equations)
    precipitation_normal = precipitation * k_normal

    return SVector(precipitation_normal,
                   precipitation_normal * (u[2] / rho),
                   precipitation_normal * (u[3] / rho),
                   precipitation_normal * (u[4] / rho),
                   zero(precipitation_normal),
                   (moisture * k_normal)...)
end

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, aux, normal_direction::AbstractVector,
                      equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta_v = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
               v3 * normal_direction[3]
    p = pressure(u, equations)

    f1 = rho * v_normal
    f2 = rho_v1 * v_normal + p * normal_direction[1]
    f3 = rho_v2 * v_normal + p * normal_direction[2]
    f4 = rho_v3 * v_normal + p * normal_direction[3]
    f5 = rho_theta_v * v_normal
    f_moisture = moisture_variables(u, equations) * v_normal

    return (SVector(f1, f2, f3, f4, f5, f_moisture...) -
            flux_sedimentation(u, aux, normal_direction, equations))
end

# Kinetic energy preserving flux of Kennedy and Gruber (2008) with the virtual potential
# temperature and the moisture variables advected as specific quantities, central fluxes of
# the precipitation, and the gravity of `flux_nonconservative_souza_etal` acting on the
# total density
@inline function flux_kennedy_gruber_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                                normal_direction::AbstractVector,
                                                equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_ll, rho_v1_ll, rho_v2_ll, rho_v3_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_v3_rr = u_rr
    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    v1_rr = rho_v1_rr / rho_rr
    v2_rr = rho_v2_rr / rho_rr
    v3_rr = rho_v3_rr / rho_rr
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2] +
                  v3_avg * normal_direction[3]

    # Conservative part
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * 0.5f0 * (u_ll[5] / rho_ll + u_rr[5] / rho_rr)
    f_moisture = f1 * 0.5f0 *
                 (moisture_variables(u_ll, equations) / rho_ll +
                  moisture_variables(u_rr, equations) / rho_rr)
    flux = (SVector(f1, f2, f3, f4, f5, f_moisture...) -
            0.5f0 * (flux_sedimentation(u_ll, aux_ll, normal_direction, equations) +
             flux_sedimentation(u_rr, aux_rr, normal_direction, equations)))

    # Nonconservative part, scaled by 0.5
    gravity = 0.5f0 * rho_avg * (aux_rr[1] - aux_ll[1])
    zero_ = zero(gravity)
    g = SVector(zero_,
                normal_direction[1] * gravity,
                normal_direction[2] * gravity,
                normal_direction[3] * gravity,
                zero_,
                ntuple(_ -> zero_, Val(nvariables(equations) - 5))...)

    return flux + g, flux - g
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_kennedy_gruber_souza_etal), ::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# `flux_ec` of the dry equations with the total density and the virtual potential
# temperature, the moisture variables advected as specific quantities with the arithmetic
# mean, central fluxes of the precipitation, and the gravity of
# `flux_nonconservative_isothermal` acting on the logarithmic mean of the total density
@inline function flux_ec(u_ll, u_rr, aux_ll, aux_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_ll, rho_v1_ll, rho_v2_ll, rho_v3_ll, rho_theta_v_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_v3_rr, rho_theta_v_rr = u_rr
    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    v1_rr = rho_v1_rr / rho_rr
    v2_rr = rho_v2_rr / rho_rr
    v3_rr = rho_v3_rr / rho_rr
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
                 v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
                 v3_rr * normal_direction[3]

    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Conservative part
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * v3_avg + p_avg * normal_direction[3]
    f5 = f1 * inv_ln_mean(rho_ll / rho_theta_v_ll, rho_rr / rho_theta_v_rr)
    f_moisture = f1 * 0.5f0 *
                 (moisture_variables(u_ll, equations) / rho_ll +
                  moisture_variables(u_rr, equations) / rho_rr)
    flux = (SVector(f1, f2, f3, f4, f5, f_moisture...) -
            0.5f0 * (flux_sedimentation(u_ll, aux_ll, normal_direction, equations) +
             flux_sedimentation(u_rr, aux_rr, normal_direction, equations)))

    # Nonconservative part, scaled by 0.5
    gravity = 0.5f0 * rho_mean * (aux_rr[1] - aux_ll[1])
    zero_ = zero(gravity)
    g = SVector(zero_,
                normal_direction[1] * gravity,
                normal_direction[2] * gravity,
                normal_direction[3] * gravity,
                zero_,
                ntuple(_ -> zero_, Val(nvariables(equations) - 5))...)

    return flux + g, flux - g
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_ec), ::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# Low Mach number approximate Riemann solver (LMARS) with the moisture variables upwinded
# with the interface velocity and without the precipitation, see `flux_lmars_combined`
@inline function flux_lmars_without_sedimentation(u_ll, u_rr, aux_ll, aux_rr,
                                                  normal_direction::AbstractVector,
                                                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

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

    p_interface = 0.5f0 * (p_ll + p_rr) -
                  0.5f0 * c_avg * rho_avg * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) -
                  1 / (2 * c_avg * rho_avg) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f = u_ll * v_interface
    else
        f = u_rr * v_interface
    end

    return SVector(f[1],
                   f[2] + p_interface * normal_direction[1],
                   f[3] + p_interface * normal_direction[2],
                   f[4] + p_interface * normal_direction[3],
                   f[5], moisture_variables(f, equations)...)
end

# LMARS with the sedimentation upwinded in the direction of the precipitation. The
# geopotential is continuous across interfaces, so the nonconservative gravity term
# vanishes on the surface and both contributions are given by the conservative flux.
@inline function flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    flux = flux_lmars_without_sedimentation(u_ll, u_rr, aux_ll, aux_rr,
                                            normal_direction, equations)

    # The precipitation falls in the direction of `-k`, i.e., from the right state if the
    # normal direction points upwards
    k_normal = 0.5f0 * (vertical_component(aux_ll, normal_direction, equations) +
                vertical_component(aux_rr, normal_direction, equations))
    if k_normal > 0
        flux = flux - flux_sedimentation(u_rr, aux_rr, normal_direction, equations)
    else
        flux = flux - flux_sedimentation(u_ll, aux_ll, normal_direction, equations)
    end

    return flux, flux
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_lmars_combined), ::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# Entropy-stable LMARS of `flux_lmars_es` of the dry equations with the total density and the
# virtual potential temperature, and without the precipitation. The moisture variables use
# the same upwinding of the mixing ratios `q = rho_x / rho` around their arithmetic mean
# instead of `1 / {{1/q}}_log`, which is not defined for `q = 0`. For `q >= 0`, both keep
# the interface value within the range of the left and right values.
@inline function flux_lmars_es_without_sedimentation(u_ll, u_rr, aux_ll, aux_rr,
                                                     normal_direction::AbstractVector,
                                                     equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

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

    # Upwinded mass and momentum fluxes
    if v_interface > 0
        f1 = u_ll[1] * v_interface
        f2 = u_ll[2] * v_interface
        f3 = u_ll[3] * v_interface
        f4 = u_ll[4] * v_interface
        sign_v = one(f1)
    else
        f1 = u_rr[1] * v_interface
        f2 = u_rr[2] * v_interface
        f3 = u_rr[3] * v_interface
        f4 = u_rr[4] * v_interface
        sign_v = -one(f1)
    end

    # Entropy-stable flux of the virtual potential temperature as in `flux_lmars_es` of the
    # dry equations
    theta_v_ll = u_ll[5] / rho_ll
    theta_v_rr = u_rr[5] / rho_rr
    theta_v_mean = inv_ln_mean(rho_ll / u_ll[5], rho_rr / u_rr[5])
    theta_v_min = min(theta_v_ll, theta_v_rr)
    theta_v_avg = 0.5f0 * (theta_v_ll + theta_v_rr)
    f5 = f1 *
         (theta_v_mean -
          theta_v_min / (2 * theta_v_avg) * (theta_v_rr - theta_v_ll) * sign_v)

    q_ll = moisture_variables(u_ll, equations) / rho_ll
    q_rr = moisture_variables(u_rr, equations) / rho_rr
    q_min = min.(q_ll, q_rr)
    # The upwinding vanishes for `q_min = 0` and is switched off for `q_min < 0`
    scaling = ifelse.(q_min .> 0, q_min ./ (q_ll + q_rr), zero(q_min))
    f_moisture = f1 * (0.5f0 * (q_ll + q_rr) - scaling .* (q_rr - q_ll) * sign_v)

    return SVector(f1,
                   f2 + p_interface * normal_direction[1],
                   f3 + p_interface * normal_direction[2],
                   f4 + p_interface * normal_direction[3],
                   f5, f_moisture...)
end

# Entropy-stable LMARS with the sedimentation upwinded as in `flux_lmars_combined`
@inline function flux_lmars_es(u_ll, u_rr, aux_ll, aux_rr,
                               normal_direction::AbstractVector,
                               equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    flux = flux_lmars_es_without_sedimentation(u_ll, u_rr, aux_ll, aux_rr,
                                               normal_direction, equations)

    # The precipitation falls in the direction of `-k`, i.e., from the right state if the
    # normal direction points upwards
    k_normal = 0.5f0 * (vertical_component(aux_ll, normal_direction, equations) +
                vertical_component(aux_rr, normal_direction, equations))
    if k_normal > 0
        flux = flux - flux_sedimentation(u_rr, aux_rr, normal_direction, equations)
    else
        flux = flux - flux_sedimentation(u_ll, aux_ll, normal_direction, equations)
    end

    return flux, flux
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_lmars_es), ::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# Surface fluxes without the precipitation, to which the slip-wall boundary condition adds
# the outflow through the lower boundary
@inline surface_flux_without_sedimentation(::typeof(flux_lmars_combined)) = flux_lmars_without_sedimentation
@inline surface_flux_without_sedimentation(::typeof(flux_lmars_es)) = flux_lmars_es_without_sedimentation

@inline function prim2cons(prim, aux,
                           equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, v1, v2, v3, p = prim
    # The primitive moisture variables are ratios with respect to the dry density
    q = moisture_variables(prim, equations)
    rho_dry = rho / (1 + density_water(q, equations.microphysics))
    rho_theta_v = (p / equations.K)^(1 / equations.gamma)
    return SVector(rho, rho * v1, rho * v2, rho * v3, rho_theta_v, (rho_dry * q)...)
end

@inline function cons2prim(u, aux,
                           equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3 = u
    rho_dry = density_dry(u, equations)
    p = pressure(u, equations)
    return SVector(rho, rho_v1 / rho, rho_v2 / rho, rho_v3 / rho, p,
                   (moisture_variables(u, equations) / rho_dry)...)
end

@inline function density_dry(u,
                             equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return u[1] -
           density_water(moisture_variables(u, equations), equations.microphysics)
end

@inline function pressure(u,
                          equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    # Dalton's law p = (rho_d R + rho_v R_v) T with T = theta * (p / p_0)^(R / c_p), i.e.,
    # p = K * (rho * theta_v)^gamma
    return equations.K * exp(equations.gamma * log(u[5]))
end

@inline function temperature(u,
                             equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_dry = density_dry(u, equations)
    return pressure(u, equations) / (rho_dry * equations.R + u[6] * equations.R_v)
end

# Largest terminal velocity of the precipitation of the state `u`
@inline function max_terminal_velocity(u,
                                       equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return max_terminal_velocity(moisture_variables(u, equations),
                                 density_dry(u, equations), equations.microphysics)
end

# The precipitation falls with the terminal velocity in the direction of `-k`
@inline function max_abs_speeds(u, aux,
                                equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, v1, v2, v3, p = cons2prim(u, aux, equations)
    c = sqrt(equations.gamma * p / rho)
    V = max_terminal_velocity(u, equations)

    return (abs(v1) + c + V * abs(aux[2]), abs(v2) + c + V * abs(aux[3]),
            abs(v3) + c + V * abs(aux[4]))
end

@inline function max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Calculate normal velocities and sound speeds
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)
    V = max(max_terminal_velocity(u_ll, equations),
            max_terminal_velocity(u_rr, equations))
    k_normal = max(abs(vertical_component(aux_ll, normal_direction, equations)),
                   abs(vertical_component(aux_rr, normal_direction, equations)))

    return (max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction) +
            V * k_normal)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, aux_ll, aux_rr,
                               normal_direction::AbstractVector,
                               equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, aux_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, aux_rr, equations)

    # Calculate normal velocities and sound speeds
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_ +
               max_terminal_velocity(u_ll, equations) *
               abs(vertical_component(aux_ll, normal_direction, equations)),
               abs(v_rr) + c_rr * norm_ +
               max_terminal_velocity(u_rr, equations) *
               abs(vertical_component(aux_rr, normal_direction, equations)))
end

@doc raw"""
    source_terms_microphysics(u, aux, x, t,
                              equations::CompressibleEulerMoistPotentialTemperatureEquations3D)

Source terms of the microphysics scheme `equations.microphysics`, e.g.,
[`KesslerMicrophysics`](@ref). The microphysics gives the right-hand sides of the moisture
variables and the isobaric heating rate ``\partial_t T``. Since the microphysics does not
change the dry air, the right-hand side of
``\varrho \theta_v = \varrho_d \theta (1 + q_v R_v / R)`` follows from
``\partial_t \theta = \partial_t T / \pi`` and the one of the vapor.
"""
@inline function source_terms_microphysics(u, aux, x, t,
                                           equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_dry = density_dry(u, equations)
    p = pressure(u, equations)
    exner = (p / equations.p_0)^(equations.R / equations.c_p)
    # Dry potential temperature theta = rho * theta_v / (rho_d * (1 + q_v R_v / R))
    virtual_factor = 1 + equations.R_v / equations.R * u[6] / rho_dry
    theta = u[5] / (rho_dry * virtual_factor)

    heating, source_moisture = microphysics_tendencies(moisture_variables(u, equations),
                                                       rho_dry, exner * theta, p,
                                                       equations.c_p,
                                                       equations.microphysics)
    source_rho_theta_v = rho_dry * virtual_factor * heating / exner +
                         equations.R_v / equations.R * theta * source_moisture[1]

    zero_ = zero(eltype(u))
    return SVector(zero_, zero_, zero_, zero_, source_rho_theta_v, source_moisture...)
end
end # @muladd

# Specializations of the `FluxTurbo` interface of Trixi.jl for the volume fluxes
# `FluxTurbo(flux_kennedy_gruber_souza_etal)` and `FluxTurbo(flux_ec)` with the
# `KesslerMicrophysics`.
# We do not wrap this code in `@muladd begin ... end` block. Optimizations like
# this are handled automatically by LoopVectorization.jl.
# We precompute `rho, v1, v2, v3, p`, the specific quantities `u[5:8] / rho`, the rain flux
# `rho_rain * V_r * k`, and the geopotential, i.e., 13 quantities. The `@turbo` kernels of
# Trixi.jl pass `2 * nturbovars + 5` arguments to `flux_turbo`, and LoopVectorization.jl
# allows at most 32, so `V_r` and `k` are not precomputed separately. Both fluxes use the
# same precomputed variables.
@inline function Trixi.nturbovars(::Union{typeof(flux_kennedy_gruber_souza_etal),
                                          typeof(flux_ec)},
                                  ::CompressibleEulerMoistPotentialTemperatureEquations3D{<:Any,
                                                                                          <:Any,
                                                                                          <:KesslerMicrophysics})
    return Val(13)
end

@inline function Trixi.cons2turbo(::Union{typeof(flux_kennedy_gruber_souza_etal),
                                          typeof(flux_ec)},
                                  rho, rho_v1, rho_v2, rho_v3, rho_theta_v, rho_vapor,
                                  rho_cloud, rho_rain, phi, k1, k2, k3,
                                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D{<:Any,
                                                                                                   <:Any,
                                                                                                   <:KesslerMicrophysics})
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = equations.K * exp(equations.gamma * log(rho_theta_v))
    # Rain flux `rho_rain * V_r` of `sedimentation`, written with `exp` and `log` for
    # LoopVectorization.jl
    RealT = typeof(rho)
    rho_dry = rho - rho_vapor - rho_cloud - rho_rain
    rain_positive = max(rho_rain, zero(RealT))
    rain = equations.microphysics.rain_flux_coefficient *
           exp(convert(RealT, 1.1364) * log(rain_positive)) *
           sqrt(equations.microphysics.rho_surface / rho_dry)

    return rho, v1, v2, v3, p, rho_theta_v / rho, rho_vapor / rho, rho_cloud / rho,
           rho_rain / rho, rain * k1, rain * k2, rain * k3, phi
end

@inline function Trixi.flux_turbo(::typeof(flux_kennedy_gruber_souza_etal),
                                  rho_ll, v1_ll, v2_ll, v3_ll, p_ll, theta_v_ll, q_v_ll,
                                  q_c_ll, q_r_ll, rain_k1_ll, rain_k2_ll, rain_k3_ll,
                                  phi_ll,
                                  rho_rr, v1_rr, v2_rr, v3_rr, p_rr, theta_v_rr, q_v_rr,
                                  q_c_rr, q_r_rr, rain_k1_rr, rain_k2_rr, rain_k3_rr,
                                  phi_rr,
                                  normal_direction_1, normal_direction_2,
                                  normal_direction_3,
                                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D{<:Any,
                                                                                                   <:Any,
                                                                                                   <:KesslerMicrophysics})
    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    v_dot_n_avg = v1_avg * normal_direction_1 + v2_avg * normal_direction_2 +
                  v3_avg * normal_direction_3

    # Central flux of the falling rain, see `flux_sedimentation`
    sedimentation_ll = rain_k1_ll * normal_direction_1 + rain_k2_ll * normal_direction_2 +
                       rain_k3_ll * normal_direction_3
    sedimentation_rr = rain_k1_rr * normal_direction_1 + rain_k2_rr * normal_direction_2 +
                       rain_k3_rr * normal_direction_3
    sedimentation_avg = 0.5f0 * (sedimentation_ll + sedimentation_rr)

    # Conservative part, i.e., `flux_kennedy_gruber_souza_etal`
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction_1 -
         0.5f0 * (sedimentation_ll * v1_ll + sedimentation_rr * v1_rr)
    f3 = f1 * v2_avg + p_avg * normal_direction_2 -
         0.5f0 * (sedimentation_ll * v2_ll + sedimentation_rr * v2_rr)
    f4 = f1 * v3_avg + p_avg * normal_direction_3 -
         0.5f0 * (sedimentation_ll * v3_ll + sedimentation_rr * v3_rr)
    f5 = f1 * 0.5f0 * (theta_v_ll + theta_v_rr)
    f6 = f1 * 0.5f0 * (q_v_ll + q_v_rr)
    f7 = f1 * 0.5f0 * (q_c_ll + q_c_rr)
    f8 = f1 * 0.5f0 * (q_r_ll + q_r_rr) - sedimentation_avg
    f1 = f1 - sedimentation_avg

    # Nonconservative part, scaled by 0.5. Since it is antisymmetric, the right contribution
    # is obtained by flipping the sign.
    gravity = 0.5f0 * rho_avg * (phi_rr - phi_ll)
    g2 = normal_direction_1 * gravity
    g3 = normal_direction_2 * gravity
    g4 = normal_direction_3 * gravity

    flux_left = SVector(f1, f2 + g2, f3 + g3, f4 + g4, f5, f6, f7, f8)
    flux_right = SVector(f1, f2 - g2, f3 - g3, f4 - g4, f5, f6, f7, f8)

    return flux_left, flux_right
end

@inline function Trixi.flux_turbo(::typeof(flux_ec),
                                  rho_ll, v1_ll, v2_ll, v3_ll, p_ll, theta_v_ll, q_v_ll,
                                  q_c_ll, q_r_ll, rain_k1_ll, rain_k2_ll, rain_k3_ll,
                                  phi_ll,
                                  rho_rr, v1_rr, v2_rr, v3_rr, p_rr, theta_v_rr, q_v_rr,
                                  q_c_rr, q_r_rr, rain_k1_rr, rain_k2_rr, rain_k3_rr,
                                  phi_rr,
                                  normal_direction_1, normal_direction_2,
                                  normal_direction_3,
                                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D{<:Any,
                                                                                                   <:Any,
                                                                                                   <:KesslerMicrophysics})
    v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
                 v3_ll * normal_direction_3
    v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
                 v3_rr * normal_direction_3

    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Central flux of the falling rain, see `flux_sedimentation`
    sedimentation_ll = rain_k1_ll * normal_direction_1 + rain_k2_ll * normal_direction_2 +
                       rain_k3_ll * normal_direction_3
    sedimentation_rr = rain_k1_rr * normal_direction_1 + rain_k2_rr * normal_direction_2 +
                       rain_k3_rr * normal_direction_3
    sedimentation_avg = 0.5f0 * (sedimentation_ll + sedimentation_rr)

    # Conservative part, i.e., `flux_ec`. Since `rho / rho_theta_v = 1 / theta_v`, the
    # virtual potential temperature flux uses the same mean value as in `flux_ec`.
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction_1 -
         0.5f0 * (sedimentation_ll * v1_ll + sedimentation_rr * v1_rr)
    f3 = f1 * v2_avg + p_avg * normal_direction_2 -
         0.5f0 * (sedimentation_ll * v2_ll + sedimentation_rr * v2_rr)
    f4 = f1 * v3_avg + p_avg * normal_direction_3 -
         0.5f0 * (sedimentation_ll * v3_ll + sedimentation_rr * v3_rr)
    f5 = f1 * inv_ln_mean(1 / theta_v_ll, 1 / theta_v_rr)
    f6 = f1 * 0.5f0 * (q_v_ll + q_v_rr)
    f7 = f1 * 0.5f0 * (q_c_ll + q_c_rr)
    f8 = f1 * 0.5f0 * (q_r_ll + q_r_rr) - sedimentation_avg
    f1 = f1 - sedimentation_avg

    # Nonconservative part, scaled by 0.5. Since it is antisymmetric, the right contribution
    # is obtained by flipping the sign.
    gravity = 0.5f0 * rho_mean * (phi_rr - phi_ll)
    g2 = normal_direction_1 * gravity
    g3 = normal_direction_2 * gravity
    g4 = normal_direction_3 * gravity

    flux_left = SVector(f1, f2 + g2, f3 + g3, f4 + g4, f5, f6, f7, f8)
    flux_right = SVector(f1, f2 - g2, f3 - g3, f4 - g4, f5, f6, f7, f8)

    return flux_left, flux_right
end
