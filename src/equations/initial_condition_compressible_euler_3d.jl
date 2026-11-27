# Unperturbed balanced steady-state.
# Returns primitive variables with only the velocity in longitudinal direction (rho, u, p).
# The other velocity components are zero.
@inline function basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)
    RealT = typeof(z)
    # Parameters from Table 1 in the paper
    # Corresponding names in the paper are commented
    radius_earth = convert(RealT, EARTH_RADIUS)                     # a
    half_width_parameter = 2                                        # b
    gravitational_acceleration = convert(RealT,
                                         EARTH_GRAVITATIONAL_ACCELERATION) # g
    k = 3                                                           # k
    surface_pressure = 1.0f5                                        # p₀
    gas_constant = 287                                              # R
    surface_equatorial_temperature = 310                            # T₀ᴱ
    surface_polar_temperature = 240                                 # T₀ᴾ
    lapse_rate = convert(RealT, 0.005)                              # Γ
    angular_velocity = convert(RealT, EARTH_ROTATION_RATE)          # Ω

    # Distance to the center of the Earth
    r = z + radius_earth

    # In the paper: T₀
    temperature0 = 0.5f0 * convert(RealT,
                           surface_equatorial_temperature + surface_polar_temperature)
    # In the paper: A, B, C, H
    const_a = 1 / lapse_rate
    const_b = (temperature0 - surface_polar_temperature) /
              (temperature0 * surface_polar_temperature)
    const_c = 0.5f0 * convert(RealT, k + 2) *
              (surface_equatorial_temperature - surface_polar_temperature) /
              (surface_equatorial_temperature * surface_polar_temperature)
    const_h = gas_constant * temperature0 / gravitational_acceleration

    # In the paper: (r - a) / bH
    scaled_z = z / (half_width_parameter * const_h)

    # Temporary variables
    temp1 = exp(lapse_rate / temperature0 * z)
    temp2 = exp(-scaled_z^2)

    # In the paper: ̃τ₁, ̃τ₂
    tau1 = const_a * lapse_rate / temperature0 * temp1 +
           const_b * (1 - 2 * scaled_z^2) * temp2
    tau2 = const_c * (1 - 2 * scaled_z^2) * temp2

    # In the paper: ∫τ₁(r') dr', ∫τ₂(r') dr'
    inttau1 = const_a * (temp1 - 1) + const_b * z * temp2
    inttau2 = const_c * z * temp2

    # Temporary variables
    temp3 = r / radius_earth * cos(lat)
    temp4 = temp3^k - convert(RealT, k) / (k + 2) * temp3^(k + 2)

    # In the paper: T
    temperature = 1 / ((r / radius_earth)^2 * (tau1 - tau2 * temp4))

    # In the paper: U, u (zonal wind, first component of spherical velocity)
    big_u = gravitational_acceleration / radius_earth * k * temperature * inttau2 *
            (temp3^(k - 1) - temp3^(k + 1))
    temp5 = radius_earth * cos(lat)
    u = -angular_velocity * temp5 + sqrt(angular_velocity^2 * temp5^2 + temp5 * big_u)

    # Hydrostatic pressure
    p = surface_pressure *
        exp(-gravitational_acceleration / gas_constant * (inttau1 - inttau2 * temp4))

    # Density (via ideal gas law)
    rho = p / (gas_constant * temperature)

    return rho, u, p
end

# Perturbation as in Equations 25 and 26 of the paper (analytical derivative)
@inline function perturbation_stream_function(lon, lat, z)
    RealT = typeof(z)
    # Parameters from Table 1 in the paper
    # Corresponding names in the paper are commented
    perturbation_radius = inv(convert(RealT, 6)) # d₀ / a
    perturbed_wind_amplitude = 1      # Vₚ
    perturbation_lon = convert(RealT, pi) / 9     # Longitude of perturbation location
    perturbation_lat = 2 * convert(RealT, pi) / 9 # Latitude of perturbation location
    pertz = 15000    # Perturbation height cap

    # Great circle distance (d in the paper) divided by a (radius of the Earth)
    # because we never actually need d without dividing by a
    great_circle_distance_by_a = acos(sin(perturbation_lat) * sin(lat) +
                                      cos(perturbation_lat) * cos(lat) *
                                      cos(lon - perturbation_lon))

    # In the first case, the vertical taper function is per definition zero.
    # In the second case, the stream function is per definition zero.
    if z > pertz || great_circle_distance_by_a > perturbation_radius
        return 0, 0
    end

    # Vertical tapering of stream function
    perttaper = 1 - 3 * z^2 / pertz^2 + 2 * z^3 / pertz^3

    # sin/cos(pi * d / (2 * d_0)) in the paper
    sin_, cos_ = sincos(pi * great_circle_distance_by_a /
                        (2 * perturbation_radius))

    # Common factor for both u and v
    factor = 16 / (3 * sqrt(convert(RealT, 3))) * perturbed_wind_amplitude * perttaper *
             cos_^3 * sin_

    u_perturbation = -factor * (-sin(perturbation_lat) * cos(lat) +
                      cos(perturbation_lat) * sin(lat) * cos(lon - perturbation_lon)) /
                     sin(great_circle_distance_by_a)

    v_perturbation = factor * cos(perturbation_lat) * sin(lon - perturbation_lon) /
                     sin(great_circle_distance_by_a)

    return u_perturbation, v_perturbation
end

@doc raw"""
    initial_condition_baroclinic_instability(x, t,
                                             equations::AbstractCompressibleEulerEquations{3, 5})

Initial condition for an idealized baroclinic instability test on a spherical shell of inner
radius [`EARTH_RADIUS`](@ref). It consists of a balanced, hydrostatic and geostrophically balanced
steady state, to which a localized stream function perturbation is added in order to trigger
the instability. The corresponding geopotential is given by
[`aux_field_gravity_earth`](@ref).

This initial condition is meant to be used together with [`source_terms_coriolis`](@ref),
which accounts for the rotation of the Earth.

See Section 3.2 and Appendix A of the following paper:
- Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013). A
  proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores.
  [DOI: 10.1002/qj.2241](https://doi.org/10.1002/qj.2241)
"""
function initial_condition_baroclinic_instability(x, t,
                                                  equations::AbstractCompressibleEulerEquations{3,
                                                                                                5})
    RealT = eltype(x)
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = convert(RealT, EARTH_RADIUS)
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Stream function type perturbation
    u_perturbation, v_perturbation = perturbation_stream_function(lon, lat, z)

    u += u_perturbation
    v = v_perturbation

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u - sin(lat) * cos(lon) * v
    v2 = cos(lon) * u - sin(lat) * sin(lon) * v
    v3 = cos(lat) * v

    aux = aux_field_gravity_earth(x, equations)
    return prim2cons(SVector(rho, v1, v2, v3, p), aux, equations)
end

@doc raw"""
    aux_field_gravity_earth(x, equations::AbstractCompressibleEulerEquations{3, 5})

Geopotential
```math
\phi(\vec{x}) = g \frac{a z}{\lVert \vec{x} \rVert}, \qquad z = \lVert \vec{x} \rVert - a,
```
of the gravity of the Earth, which points to the center of the Earth and decreases with
``1 / \lVert \vec{x} \rVert^2``. Here, ``a`` is the radius of the Earth
([`EARTH_RADIUS`](@ref)) and ``g`` the gravitational acceleration at its surface
([`EARTH_GRAVITATIONAL_ACCELERATION`](@ref)), where ``\phi = 0``, and ``\phi \approx g z``
near it. Computing it from the height ``z`` avoids the cancellation of two large terms,
which matters in `Float32`. It is returned as the auxiliary variables, see
`Trixi.nauxvariables`, and passed as `aux_field` to `SemidiscretizationHyperbolic`.
"""
@inline function aux_field_gravity_earth(x,
                                         equations::AbstractCompressibleEulerEquations{3,
                                                                                       5})
    RealT = eltype(x)
    radius_earth = convert(RealT, EARTH_RADIUS)
    gravitational_acceleration = convert(RealT, EARTH_GRAVITATIONAL_ACCELERATION)
    r = norm(x)
    phi = gravitational_acceleration * radius_earth * (r - radius_earth) / r

    return SVector(phi)
end

# Exponential perturbation of the zonal wind of the moist baroclinic wave of DCMIP2016
# with an amplitude of 1 m/s. As in DCMIP2016, it is cut off at the perturbation radius,
# where it jumps from `exp(-1)` to zero.
@inline function perturbation_exponential(lon, lat, z)
    RealT = typeof(z)
    perturbation_radius = convert(RealT, 0.1) # in units of the radius of the Earth
    perturbation_lon = convert(RealT, pi) / 9
    perturbation_lat = 2 * convert(RealT, pi) / 9
    pertz = 15000 # Perturbation height cap

    # The argument of `acos` can exceed one by round-off at the center of the perturbation
    great_circle_distance = acos(min(sin(perturbation_lat) * sin(lat) +
                                     cos(perturbation_lat) * cos(lat) *
                                     cos(lon - perturbation_lon), 1)) /
                            perturbation_radius

    if z > pertz || great_circle_distance >= 1
        return zero(RealT)
    end

    # Vertical tapering
    perttaper = 1 - 3 * z^2 / pertz^2 + 2 * z^3 / pertz^3

    return perttaper * exp(-great_circle_distance^2)
end

@doc raw"""
    initial_condition_moist_baroclinic_instability(x, t,
                                                   equations::CompressibleEulerMoistPotentialTemperatureEquations3D)

Initial condition of the moist baroclinic wave of DCMIP2016 (test case 161) in the
deep-atmosphere form on a spherical shell of inner radius [`EARTH_RADIUS`](@ref).
The pressure, the total density, and the zonal wind are the ones of the balanced state of
[`initial_condition_baroclinic_instability`](@ref), whose temperature is the virtual temperature ``T_v``. The specific humidity
``q_v = \varrho_v / \varrho`` is
```math
q_v = q_0 \exp\left(-\left(\frac{\varphi}{\varphi_w}\right)^4\right)
      \exp\left(-\left(\frac{(\eta - 1) p_0}{p_w}\right)^2\right)
```
for ``\eta = p / p_0 > 0.1`` and ``q_v = 10^{-12}`` above, with ``q_0 = 0.018``,
``\varphi_w = 2\pi / 9``, ``p_0 = 10^5`` Pa, and ``p_w = 34000`` Pa. The temperature follows
from the equation of state, i.e., it is ``T = T_v / (1 + 0.608 q_v)`` for the gas constants
``R = 287`` J/(kg K) and ``R_v = 461.5`` J/(kg K) of DCMIP2016. Initially, all other moisture
variables vanish. The instability
is triggered by the exponential perturbation of the zonal wind of DCMIP2016 with an
amplitude of 1 m/s, centered at the longitude ``\pi / 9`` and the latitude ``2 \pi / 9``.
The corresponding auxiliary variables are given by [`aux_field_gravity_earth`](@ref).

This initial condition is meant to be used together with [`source_terms_microphysics`](@ref) and
[`source_terms_coriolis`](@ref).

- Paul A. Ullrich et al. (2016)
  Dynamical Core Model Intercomparison Project (DCMIP2016) test case document
  [DOI: 10.5281/zenodo.1298671](https://doi.org/10.5281/zenodo.1298671)
"""
function initial_condition_moist_baroclinic_instability(x, t,
                                                        equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    RealT = eltype(x)
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = convert(RealT, EARTH_RADIUS)
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0)

    # Unperturbed basic state, whose density is the total density of the moist air
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    eta = p / 100_000
    if eta > convert(RealT, 0.1)
        q_v = convert(RealT, 0.018) * exp(-(lat / (2 * convert(RealT, pi) / 9))^4) *
              exp(-((eta - 1) * 100_000 / 34_000)^2)
    else
        q_v = convert(RealT, 1.0e-12)
    end

    u += perturbation_exponential(lon, lat, z)

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u
    v2 = cos(lon) * u
    v3 = zero(RealT)

    # `q_v` is the specific humidity `rho_vapor / rho`, while the primitive variables contain
    # the mixing ratio `rho_vapor / rho_dry`
    r_v = q_v / (1 - q_v)

    moisture = SVector(r_v,
                       ntuple(_ -> zero(RealT), Val(nmoisturevariables(equations) - 1))...)
    aux = aux_field_gravity_earth(x, equations)
    return prim2cons(SVector(rho, v1, v2, v3, p, moisture...), aux, equations)
end

@doc raw"""
    aux_field_gravity_earth(x,
                            equations::CompressibleEulerMoistPotentialTemperatureEquations3D)

Geopotential
```math
\phi(\vec{x}) = g \frac{a z}{\lVert \vec{x} \rVert}, \qquad z = \lVert \vec{x} \rVert - a,
```
of the gravity of the Earth, which points to the center of the Earth and decreases with
``1 / \lVert \vec{x} \rVert^2``, and the upward unit vector
``\hat{k} = \vec{x} / \lVert \vec{x} \rVert``. Here, ``a`` is the radius of the Earth
([`EARTH_RADIUS`](@ref)) and ``g`` the gravitational acceleration at its surface
([`EARTH_GRAVITATIONAL_ACCELERATION`](@ref)), where ``\phi = 0``, and ``\phi \approx g z``
near it. Computing it from the height ``z`` avoids the cancellation of two large terms,
which matters in `Float32`. They are returned as the auxiliary variables
`SVector(phi, k1, k2, k3)`, see `Trixi.nauxvariables`, and passed as `aux_field` to
`SemidiscretizationHyperbolic`.
"""
@inline function aux_field_gravity_earth(x,
                                         equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    RealT = eltype(x)
    radius_earth = convert(RealT, EARTH_RADIUS)
    gravitational_acceleration = convert(RealT, EARTH_GRAVITATIONAL_ACCELERATION)
    r = norm(x)
    phi = gravitational_acceleration * radius_earth * (r - radius_earth) / r

    return SVector(phi, x[1] / r, x[2] / r, x[3] / r)
end

@doc raw"""
    initial_condition_isothermal(x, t,
                                 equations::AbstractCompressibleEulerEquations{3, 5})

Isothermal atmosphere at rest in hydrostatic balance on a spherical shell of inner radius
[`EARTH_RADIUS`](@ref), used as the initial state of the Held-Suarez test case. The
corresponding geopotential is given by [`aux_field_gravity_earth`](@ref).

- Souza et al. (2023). The Flux-Differencing Discontinuous Galerkin Method Applied to an
  Idealized Fully Compressible Nonhydrostatic Dry Atmosphere.
  [DOI: 10.1029/2022MS003527](https://doi.org/10.1029/2022MS003527)
"""
function initial_condition_isothermal(x, t,
                                      equations::AbstractCompressibleEulerEquations{3, 5})
    # equation (60) in the paper
    T = 285

    @unpack p_0, R = equations

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - EARTH_RADIUS, 0)
    r = z + EARTH_RADIUS

    # pressure, geopotential formulation
    p = p_0 *
        exp(-EARTH_GRAVITATIONAL_ACCELERATION * EARTH_RADIUS * z / r / (R * T))

    # density (via ideal gas law)
    rho = p / (R * T)

    aux = aux_field_gravity_earth(x, equations)
    return prim2cons(SVector(rho, 0, 0, 0, p), aux, equations)
end

@doc raw"""
    initial_condition_vortex_shedding(x, t,
                                      equations::AbstractCompressibleEulerEquations{3, 5})

Initial condition for an idealized mountain-triggered mesoscale flow test case (DCMIP-2025
Test Case 2), which examines the interaction between a stratified rotating flow and an
isolated topography, leading to vortex shedding in the lee of the mountain. It consists of a
stably stratified atmosphere with constant Brunt-Väisälä frequency
``N = 0.0182\ \mathrm{s}^{-1}`` in solid body rotation with
``u_0 = 20\ \mathrm{m}/\mathrm{s}``.

The test case uses the reduced-radius Earth configuration described by
[`SMALL_EARTH_SCALING`](@ref), and is therefore meant to be used together with
[`source_terms_coriolis_small_earth`](@ref). The corresponding geopotential is given by
[`aux_field_gravity_small_earth`](@ref).

- DCMIP-2025 Organizing Committee (2025). Test Case 2: Mountain-Triggered Meso-scale Flow
  Phenomena.
  [https://sites.google.com/umich.edu/dcmip-2025](https://sites.google.com/umich.edu/dcmip-2025/dcmip-2025-test-cases/test-case-2-mountain-triggered-meso-scale-flow-phenomena)
"""
function initial_condition_vortex_shedding(x, t,
                                           equations::AbstractCompressibleEulerEquations{3,
                                                                                         5})
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = SMALL_EARTH_RADIUS
    R = equations.c_p - equations.c_v
    k = R / equations.c_p
    # Convert spherical velocity to Cartesian
    u0 = 20
    v1 = -u0 * cos(lat) * sin(lon)
    v2 = u0 * cos(lat) * cos(lon)
    v3 = 0
    g = EARTH_GRAVITATIONAL_ACCELERATION  # g
    T = 288
    angular_velocity = SMALL_EARTH_ROTATION_RATE  # Ω
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0)
    phi = g * z
    N = 0.0182 # Brunt–Väisälä frequency
    p = equations.p_0 * exp(-radius_earth * N^2 * u0 / (2 * g^2 * k) *
            (u0 / radius_earth + 2 * angular_velocity) * (sin(lat)^2 - 1) -
            N^2 / (g^2 * k) * phi)
    rho = p / (R * T)
    aux = aux_field_gravity_small_earth(x, equations)
    return prim2cons(SVector(rho, v1, v2, v3, p), aux, equations)
end

@doc raw"""
    aux_field_gravity_small_earth(x,
                                  equations::AbstractCompressibleEulerEquations{3, 5})

Geopotential
```math
\phi(\vec{x}) = g \left(\lVert \vec{x} \rVert - \frac{a}{X}\right)
```
of a gravity of constant magnitude ``g`` ([`EARTH_GRAVITATIONAL_ACCELERATION`](@ref))
pointing to the center of the reduced-radius Earth of radius ``a / X``
([`SMALL_EARTH_RADIUS`](@ref)), see [`SMALL_EARTH_SCALING`](@ref). It is returned as the
auxiliary variables, see `Trixi.nauxvariables`, and passed as `aux_field` to
`SemidiscretizationHyperbolic`.
"""
@inline function aux_field_gravity_small_earth(x,
                                               equations::AbstractCompressibleEulerEquations{3,
                                                                                             5})
    RealT = eltype(x)
    radius_earth = convert(RealT, SMALL_EARTH_RADIUS)
    gravitational_acceleration = convert(RealT, EARTH_GRAVITATIONAL_ACCELERATION)
    phi = gravitational_acceleration * (norm(x) - radius_earth)

    return SVector(phi)
end

# Coriolis force for a rotation about the z axis with the given angular velocity, used by
# `source_terms_coriolis` and `source_terms_coriolis_small_earth`
@inline function coriolis_source(u, aux, angular_velocity,
                                 equations::AbstractCompressibleEulerEquations{3})
    du0 = zero(eltype(u))
    angular_velocity = convert(eltype(u), angular_velocity)

    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 = 2 * angular_velocity * u[3]
    du3 = -2 * angular_velocity * u[2]

    return SVector(du0, du2, du3, du0,
                   ntuple(_ -> du0, Val(nvariables(equations) - 4))...)
end

@doc raw"""
    source_terms_coriolis(u, aux, x, t,
                          equations::AbstractCompressibleEulerEquations{3, 5})

Source term accounting for the Coriolis force ``-2 \vec{\Omega} \times \rho \vec{v}``
induced by the rotation of the Earth about the ``z`` axis, with the rotation rate given by
[`EARTH_ROTATION_RATE`](@ref). For the [`CompressibleEulerVectorInvariantEquations3D`](@ref),
whose variables contain the velocity instead of the momentum, this is the Coriolis force
``-2 \vec{\Omega} \times \vec{v}``.

See also [`source_terms_coriolis_small_earth`](@ref) for the reduced-radius configuration.
"""
@inline function source_terms_coriolis(u, aux, x, t,
                                       equations::AbstractCompressibleEulerEquations{3, 5})
    return coriolis_source(u, aux, EARTH_ROTATION_RATE, equations)
end

# The Coriolis force acts on the momentum of the total density
@inline function source_terms_coriolis(u, aux, x, t,
                                       equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return coriolis_source(u, aux, EARTH_ROTATION_RATE, equations)
end

@doc raw"""
    source_terms_coriolis_small_earth(u, aux, x, t,
                                      equations::AbstractCompressibleEulerEquations{3, 5})

Counterpart of [`source_terms_coriolis`](@ref) for the reduced-radius ("small Earth")
configuration, using the scaled rotation rate [`SMALL_EARTH_ROTATION_RATE`](@ref).
"""
@inline function source_terms_coriolis_small_earth(u, aux, x, t,
                                                   equations::AbstractCompressibleEulerEquations{3,
                                                                                                 5})
    return coriolis_source(u, aux, SMALL_EARTH_ROTATION_RATE, equations)
end
