@muladd begin
#! format: noindent

@doc raw"""
    CompressibleEulerMoistPotentialTemperatureEquations3D(; c_p, c_v, gravity,
                                                          R_v = 461.5,
                                                          rho_surface = 1.15,
                                                          condensation_time = 1)

The compressible Euler equations with gravity in potential temperature formulation with the
warm-rain moisture of the Kessler microphysics, as used in the test cases of DCMIP2016,
i.e., water vapor, cloud water, and rain water,
```math
\begin{aligned}
\partial_t \varrho + \nabla \cdot (\varrho \vec{v} - \varrho_r V_r \hat{k}) &= 0, \\
\partial_t (\varrho \vec{v}) + \nabla \cdot (\varrho \vec{v} \otimes \vec{v} + p I
  - \varrho_r V_r \vec{v} \otimes \hat{k}) &= - \varrho \nabla \phi, \\
\partial_t (\varrho \theta_v) + \nabla \cdot (\varrho \theta_v \vec{v})
  &= \varrho_d \left(1 + \frac{R_v}{R} q_v\right) \frac{L}{c_p \pi} (C - E)
  + \varrho_d \frac{R_v}{R} \theta (E - C), \\
\partial_t \varrho_v + \nabla \cdot (\varrho_v \vec{v}) &= \varrho_d (E - C), \\
\partial_t \varrho_c + \nabla \cdot (\varrho_c \vec{v}) &= \varrho_d (C - A - C_r), \\
\partial_t \varrho_r + \nabla \cdot (\varrho_r \vec{v} - \varrho_r V_r \hat{k})
  &= \varrho_d (A + C_r - E),
\end{aligned}
```
for the total density ``\varrho = \varrho_d + \varrho_v + \varrho_c + \varrho_r``, the
densities of dry air, vapor, cloud water, and rain water, and the virtual potential
temperature
```math
\theta_v = \frac{p}{\varrho R \pi} = \theta \frac{1 + q_v R_v / R}{1 + q_v + q_c + q_r},
```
which includes the weight of the cloud and rain water, with the dry potential temperature
``\theta = T / \pi``, the Exner pressure ``\pi = (p / p_0)^{R / c_p}``, and the mixing
ratios ``q_x = \varrho_x / \varrho_d``. The rain falls in the direction ``-\hat{k}`` with the terminal velocity ``V_r`` of
Klemp and Wilhelmson (1978) (see [`terminal_velocity_rain`](@ref)), where ``\hat{k}`` is
the upward unit vector and the density `rho_surface` replaces the dry density at the
surface. The pressure follows from Dalton's law ``p = (\varrho_d R + \varrho_v R_v) T``,
i.e., ``p = K (\varrho \theta_v)^\gamma``, which depends on the fifth variable only.

The auxiliary variables, see `Trixi.nauxvariables`, are the geopotential ``\phi`` and the
Cartesian components of ``\hat{k}``, i.e., `aux = SVector(phi, k1, k2, k3)`. Their values
are given by the function `aux_field(x, equations)` passed to
`SemidiscretizationHyperbolic`. On the sphere, ``\hat{k} = \vec{x} / \lVert \vec{x} \rVert``
is the radial direction, which varies within each element. The pointwise functions take the
auxiliary variables as additional argument after the state, e.g.,
`flux(u, aux, normal_direction, equations)`. The numerical fluxes return the conservative
and the nonconservative contribution together, see
[`Trixi.combine_conservative_and_nonconservative_fluxes`](@ref), i.e., use
[`flux_kennedy_gruber_souza_etal`](@ref) as `volume_flux` and
[`flux_lmars_combined`](@ref) as `surface_flux`.

The right-hand sides are the warm-rain microphysics of Kessler (1969) in the form of Klemp
and Wilhelmson (1978) as in DCMIP2016, with the condensation `C`, the autoconversion `A`,
the accretion `C_r`, and the evaporation of rain `E`, see [`source_terms_kessler`](@ref).
The condensation relaxes the vapor to saturation with the time scale `condensation_time`
in s.

## References

- Edwin Kessler (1969)
  On the distribution and continuity of water substance in atmospheric circulations
  [DOI: 10.1007/978-1-935704-36-2_1](https://doi.org/10.1007/978-1-935704-36-2_1)
- Joseph B. Klemp, Robert B. Wilhelmson (1978)
  The simulation of three-dimensional convective storm dynamics
  [DOI: 10.1175/1520-0469(1978)035<1070:TSOTDC>2.0.CO;2](https://doi.org/10.1175/1520-0469(1978)035<1070:TSOTDC>2.0.CO;2)
- Paul A. Ullrich et al. (2016)
  Dynamical Core Model Intercomparison Project (DCMIP2016) test case document
  [DOI: 10.5281/zenodo.1298671](https://doi.org/10.5281/zenodo.1298671)
"""
struct CompressibleEulerMoistPotentialTemperatureEquations3D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{3, 8}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat of dry air at constant pressure in J/(kg K)
    c_v::RealT # specific heat of dry air at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant of dry air in J/(kg K)
    R_v::RealT # gas constant of water vapor in J/(kg K)
    gamma::RealT # ratio of specific heats
    # Scaling factor between the pressure and the weighted potential temperature
    K::RealT # = p_0 * (R / p_0)^gamma
    rho_surface::RealT # density at the surface in the terminal velocity of rain in kg/m³
    condensation_time::RealT # relaxation time of the condensation in s
    function CompressibleEulerMoistPotentialTemperatureEquations3D(; c_p, c_v, gravity,
                                                                   R_v = 461.5,
                                                                   rho_surface = 1.15,
                                                                   condensation_time = 1)
        c_p, c_v, g = promote(c_p, c_v, gravity)
        R_v, rho_surface, condensation_time = convert.(typeof(c_p),
                                                       (R_v, rho_surface,
                                                        condensation_time))
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        K = p_0 * (R / p_0)^gamma
        return new{typeof(c_p)}(p_0, c_p, c_v, g, R, R_v, gamma, K, rho_surface,
                                condensation_time)
    end
end

# Together with the specialization of `Adapt.adapt_structure` in Trixi.jl,
# this allows to move semidiscretizations and their components including
# the equations to GPUs and adapt the floating point type, e.g.,
# to `Float32` to improve performance on GPUs.
function Base.similar(equations::CompressibleEulerMoistPotentialTemperatureEquations3D,
                      ::Type{NewRealT}) where {NewRealT}
    c_p, c_v, gravity, R_v, rho_surface, condensation_time = convert.(NewRealT,
                                                                      (equations.c_p,
                                                                       equations.c_v,
                                                                       equations.g,
                                                                       equations.R_v,
                                                                       equations.rho_surface,
                                                                       equations.condensation_time))
    return CompressibleEulerMoistPotentialTemperatureEquations3D(; c_p, c_v, gravity,
                                                                 R_v,
                                                                 rho_surface,
                                                                 condensation_time)
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerMoistPotentialTemperatureEquations3D)
    ("rho", "rho_v1", "rho_v2", "rho_v3", "rho_theta_v", "rho_vapor", "rho_cloud",
     "rho_rain")
end

# The primitive variables use the mixing ratios with respect to dry air, as in DCMIP2016
function varnames(::typeof(cons2prim),
                  ::CompressibleEulerMoistPotentialTemperatureEquations3D)
    ("rho", "v1", "v2", "v3", "p", "q_vapor", "q_cloud", "q_rain")
end

@inline Trixi.nauxvariables(::CompressibleEulerMoistPotentialTemperatureEquations3D) = 4

Trixi.auxvarnames(::CompressibleEulerMoistPotentialTemperatureEquations3D) = ("phi",
                                                                              "k1",
                                                                              "k2",
                                                                              "k3")

have_nonconservative_terms(::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# Slip-wall boundary condition
# As usual, the external state mirrors the normal velocity. In addition, it contains no
# rain: `flux_lmars_combined` upwinds the falling rain, i.e., it takes the rain of the state
# above the face. At the lower boundary, this is the inner state, so that the rain leaves
# the domain (precipitation), and at the upper boundary, it is the external state, so that
# no rain enters the domain. The rain is removed from the total density, such that the
# external state has the same dry air, pressure, and velocity as the inner one. There are
# no separate auxiliary variables outside of the domain, so the inner ones are used on both
# sides.
@inline function boundary_condition_slip_wall(u_inner, aux_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_function::typeof(flux_lmars_combined),
                                              equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    rho, rho_v1, rho_v2, rho_v3, rho_theta_v, rho_vapor, rho_cloud, rho_rain = u_inner
    u_normal = normal[1] * rho_v1 + normal[2] * rho_v2 + normal[3] * rho_v3

    # create the "external" boundary solution state without rain
    rho_boundary = rho - rho_rain
    scaling = rho_boundary / rho
    u_boundary = SVector(rho_boundary,
                         scaling * (rho_v1 - 2 * u_normal * normal[1]),
                         scaling * (rho_v2 - 2 * u_normal * normal[2]),
                         scaling * (rho_v3 - 2 * u_normal * normal[3]),
                         rho_theta_v, rho_vapor, rho_cloud, zero(rho_rain))

    # calculate the boundary flux
    flux, _ = surface_flux_function(u_inner, u_boundary, aux_inner, aux_inner,
                                    normal_direction, equations)
    return flux
end

# Component of the `normal_direction` in the direction of the upward unit vector `k`, which
# is given by the auxiliary variables `aux[2:4]` in the physical coordinates
@inline function vertical_component(aux, normal_direction::AbstractVector)
    return aux[2] * normal_direction[1] + aux[3] * normal_direction[2] +
           aux[4] * normal_direction[3]
end

# Flux of the falling rain through a face with the normal direction `normal_direction`,
# i.e., of the rain water and of its momentum, for the velocity `(v1, v2, v3)` of the
# state `u`
@inline function flux_sedimentation(u, aux, v1, v2, v3,
                                    normal_direction::AbstractVector,
                                    equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rain = u[8] * terminal_velocity_rain(u, equations) *
           vertical_component(aux, normal_direction)
    return rain, rain * v1, rain * v2, rain * v3
end

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, aux, normal_direction::AbstractVector,
                      equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, rho_theta_v, rho_vapor, rho_cloud, rho_rain = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
               v3 * normal_direction[3]
    p = pressure(u, equations)
    rain, rain_v1, rain_v2, rain_v3 = flux_sedimentation(u, aux, v1, v2, v3,
                                                         normal_direction, equations)

    f1 = rho * v_normal - rain
    f2 = rho_v1 * v_normal + p * normal_direction[1] - rain_v1
    f3 = rho_v2 * v_normal + p * normal_direction[2] - rain_v2
    f4 = rho_v3 * v_normal + p * normal_direction[3] - rain_v3
    f5 = rho_theta_v * v_normal
    f6 = rho_vapor * v_normal
    f7 = rho_cloud * v_normal
    f8 = rho_rain * v_normal - rain

    return SVector(f1, f2, f3, f4, f5, f6, f7, f8)
end

# Kinetic energy preserving flux of Kennedy and Gruber (2008) with the potential temperature
# and the water species advected as specific quantities, central fluxes of the falling rain,
# and the gravity of `flux_nonconservative_souza_etal` acting on the total density
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
    rain_ll, rain_v1_ll, rain_v2_ll, rain_v3_ll = flux_sedimentation(u_ll, aux_ll,
                                                                     v1_ll, v2_ll,
                                                                     v3_ll,
                                                                     normal_direction,
                                                                     equations)
    rain_rr, rain_v1_rr, rain_v2_rr, rain_v3_rr = flux_sedimentation(u_rr, aux_rr,
                                                                     v1_rr, v2_rr,
                                                                     v3_rr,
                                                                     normal_direction,
                                                                     equations)
    rain_avg = 0.5f0 * (rain_ll + rain_rr)

    # Conservative part
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1] - 0.5f0 * (rain_v1_ll + rain_v1_rr)
    f3 = f1 * v2_avg + p_avg * normal_direction[2] - 0.5f0 * (rain_v2_ll + rain_v2_rr)
    f4 = f1 * v3_avg + p_avg * normal_direction[3] - 0.5f0 * (rain_v3_ll + rain_v3_rr)
    f5 = f1 * 0.5f0 * (u_ll[5] / rho_ll + u_rr[5] / rho_rr)
    f6 = f1 * 0.5f0 * (u_ll[6] / rho_ll + u_rr[6] / rho_rr)
    f7 = f1 * 0.5f0 * (u_ll[7] / rho_ll + u_rr[7] / rho_rr)
    f8 = f1 * 0.5f0 * (u_ll[8] / rho_ll + u_rr[8] / rho_rr) - rain_avg
    f1 = f1 - rain_avg

    # Nonconservative part, scaled by 0.5
    gravity = 0.5f0 * rho_avg * (aux_rr[1] - aux_ll[1])
    g2 = normal_direction[1] * gravity
    g3 = normal_direction[2] * gravity
    g4 = normal_direction[3] * gravity

    return SVector(f1, f2 + g2, f3 + g3, f4 + g4, f5, f6, f7, f8),
           SVector(f1, f2 - g2, f3 - g3, f4 - g4, f5, f6, f7, f8)
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_kennedy_gruber_souza_etal), ::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

# Low Mach number approximate Riemann solver (LMARS) with the water species upwinded with
# the interface velocity and the sedimentation upwinded in the direction of the falling
# rain. The geopotential is continuous across interfaces, so the nonconservative gravity
# term vanishes on the surface and both contributions are given by the conservative flux.
@inline function flux_lmars_combined(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

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
        f1, f2, f3, f4, f5, f6, f7, f8 = u_ll * v_interface
    else
        f1, f2, f3, f4, f5, f6, f7, f8 = u_rr * v_interface
    end

    # The rain falls in the direction of `-k`, i.e., from the right state if the normal
    # direction points upwards
    k_normal = 0.5f0 * (vertical_component(aux_ll, normal_direction) +
                vertical_component(aux_rr, normal_direction))
    if k_normal > 0
        rain, rain_v1, rain_v2, rain_v3 = flux_sedimentation(u_rr, aux_rr,
                                                             v1_rr, v2_rr, v3_rr,
                                                             normal_direction,
                                                             equations)
    else
        rain, rain_v1, rain_v2, rain_v3 = flux_sedimentation(u_ll, aux_ll,
                                                             v1_ll, v2_ll, v3_ll,
                                                             normal_direction,
                                                             equations)
    end

    flux = SVector(f1 - rain,
                   f2 + p_interface * normal_direction[1] - rain_v1,
                   f3 + p_interface * normal_direction[2] - rain_v2,
                   f4 + p_interface * normal_direction[3] - rain_v3,
                   f5, f6, f7, f8 - rain)

    return flux, flux
end

@inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_lmars_combined), ::CompressibleEulerMoistPotentialTemperatureEquations3D) = Trixi.True()

@inline function prim2cons(prim,
                           equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, v1, v2, v3, p, q_v, q_c, q_r = prim
    rho_dry = rho / (1 + q_v + q_c + q_r)
    rho_theta_v = (p / equations.K)^(1 / equations.gamma)
    return SVector(rho, rho * v1, rho * v2, rho * v3, rho_theta_v,
                   rho_dry * q_v, rho_dry * q_c, rho_dry * q_r)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, rho_v1, rho_v2, rho_v3, _, rho_vapor, rho_cloud, rho_rain = u
    rho_dry = density_dry(u, equations)
    p = pressure(u, equations)
    return SVector(rho, rho_v1 / rho, rho_v2 / rho, rho_v3 / rho, p,
                   rho_vapor / rho_dry, rho_cloud / rho_dry, rho_rain / rho_dry)
end

@inline function cons2cons(u,
                           equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return u
end

@inline function density_dry(u,
                             equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, _, _, _, _, rho_vapor, rho_cloud, rho_rain = u
    return rho - rho_vapor - rho_cloud - rho_rain
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

@doc raw"""
    terminal_velocity_rain(u, equations::CompressibleEulerMoistPotentialTemperatureEquations3D)

Terminal velocity of rain in m/s of Klemp and Wilhelmson (1978, Eq. 2.15) as in DCMIP2016,
```math
V_r = 36.34 \left(10^{-3} \varrho_r\right)^{0.1364} \sqrt{\frac{\varrho_{sfc}}{\varrho_d}},
```
with the rain water density ``\varrho_r = \varrho_d q_r`` in kg/m³, i.e.,
``10^{-3} \varrho_r`` in g/cm³, and the density `rho_surface` of the equations as
``\varrho_{sfc}``.
"""
@inline function terminal_velocity_rain(u,
                                        equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    RealT = eltype(u)
    rho_rain = max(u[8], zero(RealT))
    return convert(RealT, 36.34) *
           (convert(RealT, 1.0e-3) * rho_rain)^convert(RealT, 0.1364) *
           sqrt(equations.rho_surface / density_dry(u, equations))
end

# The rain falls with the terminal velocity in the direction of `-k`
@inline function max_abs_speeds(u, aux,
                                equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho, v1, v2, v3, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)
    V_r = terminal_velocity_rain(u, equations)

    return (abs(v1) + c + V_r * abs(aux[2]), abs(v2) + c + V_r * abs(aux[3]),
            abs(v3) + c + V_r * abs(aux[4]))
end

@inline function max_abs_speed_naive(u_ll, u_rr, aux_ll, aux_rr,
                                     normal_direction::AbstractVector,
                                     equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speeds
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)
    V_r = max(terminal_velocity_rain(u_ll, equations),
              terminal_velocity_rain(u_rr, equations))
    k_normal = max(abs(vertical_component(aux_ll, normal_direction)),
                   abs(vertical_component(aux_rr, normal_direction)))

    return (max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction) +
            V_r * k_normal)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, aux_ll, aux_rr,
                               normal_direction::AbstractVector,
                               equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speeds
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] +
           v3_ll * normal_direction[3]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] +
           v3_rr * normal_direction[3]
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_ +
               terminal_velocity_rain(u_ll, equations) *
               abs(vertical_component(aux_ll, normal_direction)),
               abs(v_rr) + c_rr * norm_ +
               terminal_velocity_rain(u_rr, equations) *
               abs(vertical_component(aux_rr, normal_direction)))
end

###############################################################################
# Kessler microphysics of DCMIP2016

@doc raw"""
    saturation_mixing_ratio_tetens(T, p)

Saturation mixing ratio of water vapor of Tetens as in Klemp and Wilhelmson
(1978, Eq. 2.11) and DCMIP2016, with the temperature `T` in K and the pressure `p` in Pa,
```math
q_{vs}(T, p) = \frac{380}{p} \exp\left(17.27 \frac{T - 273}{T - 36}\right).
```
"""
@inline function saturation_mixing_ratio_tetens(T, p)
    RealT = typeof(T)
    return convert(RealT, 380) / p * exp(convert(RealT, 17.27) * (T - 273) / (T - 36))
end

# Latent heat of vaporization divided by the specific heat of dry air at constant pressure
# with the constants of DCMIP2016, i.e., `L = 2.5e6` J/kg and `c_p = 1003` J/(kg K)
@inline latent_heat_over_specific_heat(RealT) = convert(RealT, 2.5e6 / 1003)

@doc raw"""
    source_terms_kessler(u, aux, x, t,
                         equations::CompressibleEulerMoistPotentialTemperatureEquations3D)

Source terms of the warm-rain microphysics of Kessler (1969) in the form of Klemp and
Wilhelmson (1978, KW) as in DCMIP2016, with the mixing ratios
``q_x = \varrho_x / \varrho_d``,
- the condensation (``C > 0``) or evaporation of cloud water (``C < 0``),
  ```math
  C = \frac{1}{\tau} \max\left(\frac{q_v - q_{vs}}
      {1 + q_{vs} \frac{237.3 \cdot 17.27 L / c_p}{(T - 36)^2}}, -q_c\right),
  ```
  i.e., the saturation adjustment of KW (3.10) relaxed with the time scale
  ``\tau`` `= equations.condensation_time`, with the saturation mixing ratio ``q_{vs}`` of
  [`saturation_mixing_ratio_tetens`](@ref),
- the autoconversion ``A = 10^{-3} \max(q_c - 10^{-3}, 0)`` and the accretion
  ``C_r = 2.2 q_c q_r^{0.875}`` of cloud water to rain (KW 2.13),
- the evaporation of rain (KW 2.14), with the density ``r = 10^{-3} \varrho_d`` in g/cm³ and
  the pressure in mb,
  ```math
  E = \frac{(1.6 + 124.9 (r q_r)^{0.2046}) (r q_r)^{0.525}}
      {2.55 \cdot 10^6 / (p q_{vs}) + 5.4 \cdot 10^5}
      \frac{\max(q_{vs} - q_v, 0)}{r q_{vs}},
  ```
and the isobaric latent heating ``L / (c_p \pi) (C - E)`` of the dry potential temperature
with the constants ``L = 2.5 \cdot 10^6`` J/kg and ``c_p = 1003`` J/(kg K) of DCMIP2016. The
right-hand side of ``\varrho \theta_v = \varrho_d \theta (1 + q_v R_v / R)`` follows from
the ones of ``\varrho_d \theta`` and of the vapor, since the microphysics does not change the
dry air.
"""
@inline function source_terms_kessler(u, aux, x, t,
                                      equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    RealT = eltype(u)
    zero_ = zero(RealT)
    rho_vapor, rho_cloud, rho_rain = u[6], u[7], u[8]
    rho_dry = density_dry(u, equations)
    q_v = max(rho_vapor, zero_) / rho_dry
    q_c = max(rho_cloud, zero_) / rho_dry
    q_r = max(rho_rain, zero_) / rho_dry
    p = pressure(u, equations)
    exner = (p / equations.p_0)^(equations.R / equations.c_p)
    # Dry potential temperature theta = rho * theta_v / (rho_d * (1 + q_v R_v / R))
    virtual_factor = 1 + equations.R_v / equations.R * rho_vapor / rho_dry
    theta = u[5] / (rho_dry * virtual_factor)
    T = exner * theta
    latent_heat_over_c_p = latent_heat_over_specific_heat(RealT)

    # Condensation, where the denominator accounts for the increase of the saturation
    # mixing ratio due to the latent heating, with the constant `237.3` of `kessler.f90` of
    # DCMIP2016 instead of the exact value `237 = 273 - 36`
    q_vs = saturation_mixing_ratio_tetens(T, p)
    condensation = max((q_v - q_vs) /
                       (1 +
                        q_vs * convert(RealT, 237.3 * 17.27) * latent_heat_over_c_p /
                        (T - 36)^2), -q_c) / equations.condensation_time

    # Autoconversion and accretion
    autoconversion = convert(RealT, 1.0e-3) * max(q_c - convert(RealT, 1.0e-3), zero_)
    accretion = convert(RealT, 2.2) * q_c * q_r^convert(RealT, 0.875)

    # Evaporation of rain
    r = convert(RealT, 1.0e-3) * rho_dry
    rain = r * q_r
    evaporation = ((convert(RealT, 1.6) +
                    convert(RealT, 124.9) * rain^convert(RealT, 0.2046)) *
                   rain^convert(RealT, 0.525)) /
                  (convert(RealT, 2.55e6) / (convert(RealT, 1.0e-2) * p * q_vs) +
                   convert(RealT, 5.4e5)) *
                  max(q_vs - q_v, zero_) / (r * q_vs)

    source_rho_theta = rho_dry * latent_heat_over_c_p / exner *
                       (condensation - evaporation)
    source_rho_vapor = rho_dry * (evaporation - condensation)
    return SVector(zero_, zero_, zero_, zero_,
                   virtual_factor * source_rho_theta +
                   equations.R_v / equations.R * theta * source_rho_vapor,
                   source_rho_vapor,
                   rho_dry * (condensation - autoconversion - accretion),
                   rho_dry * (autoconversion + accretion - evaporation))
end
end # @muladd

# Specialization of the `FluxTurbo` interface of Trixi.jl for the volume flux
# `FluxTurbo(flux_kennedy_gruber_souza_etal)`.
# We do not wrap this code in `@muladd begin ... end` block. Optimizations like
# this are handled automatically by LoopVectorization.jl.
# We precompute `rho, v1, v2, v3, p`, the specific quantities `u[5:8] / rho`, the rain flux
# `rho_rain * V_r * k`, and the geopotential, i.e., 13 quantities, since the kernels of
# LoopVectorization.jl take at most 32 arguments.
@inline function Trixi.nturbovars(::typeof(flux_kennedy_gruber_souza_etal),
                                  ::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return Val(13)
end

@inline function Trixi.cons2turbo(::typeof(flux_kennedy_gruber_souza_etal),
                                  rho, rho_v1, rho_v2, rho_v3, rho_theta_v, rho_vapor,
                                  rho_cloud, rho_rain, phi, k1, k2, k3,
                                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = equations.K * exp(equations.gamma * log(rho_theta_v))
    # Rain flux `rho_rain * V_r` with the terminal velocity of `terminal_velocity_rain`,
    # written with `exp` and `log` for LoopVectorization.jl
    RealT = typeof(rho)
    rho_dry = rho - rho_vapor - rho_cloud - rho_rain
    rain_positive = max(rho_rain, zero(RealT))
    rain = rho_rain * convert(RealT, 36.34) *
           exp(convert(RealT, 0.1364) * log(convert(RealT, 1.0e-3) * rain_positive)) *
           sqrt(equations.rho_surface / rho_dry)

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
                                  equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
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
