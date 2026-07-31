"""
	initial_condition_gravity_waves(x, t,
                                        equations::AbstractCompressibleEulerEquations{2, 5})

Test cases for linearized analytical solution by
-  Baldauf, Michael and Brdar, Slavko (2013)
   An analytic solution for linear gravity waves in a channel as a test
   for numerical models using the non-hydrostatic, compressible {E}uler equations
   [DOI: 10.1002/qj.2105] (https://doi.org/10.1002/qj.2105)
"""
function initial_condition_gravity_waves(x, t,
                                         equations::AbstractCompressibleEulerEquations{2,
                                                                                       5})
    g = equations.g
    c_p = equations.c_p
    c_v = equations.c_v
    # center of perturbation
    x_c = 100_000
    a = 5_000
    H = 10_000
    R = c_p - c_v    # gas constant (dry air)
    T0 = 250
    delta = g / (R * T0)
    DeltaT = 0.001
    Tb = DeltaT * sinpi(x[2] / H) * exp(-(x[1] - x_c)^2 / a^2)
    ps = 100_000  # reference pressure
    rhos = ps / (T0 * R)
    rho_b = rhos * (-Tb / T0)
    p = ps * exp(-delta * x[2])
    rho = rhos * exp(-delta * x[2]) + rho_b * exp(-0.5 * delta * x[2])
    v1 = 20
    v2 = 0

    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

@doc raw"""
    initial_condition_adiabatic(x, t,
                                equations::AbstractCompressibleEulerEquations{2, 5})

Atmosphere at rest with constant potential temperature ``\theta = 300\ \mathrm{K}`` in
hydrostatic balance, used to assess the well-balanced property of a discretization. 
"""
function initial_condition_adiabatic(x, t,
                                     equations::AbstractCompressibleEulerEquations{2, 5})
    g = equations.g
    c_p = equations.c_p
    c_v = equations.c_v
    p0 = 100_000
    R = c_p - c_v    # gas constant (dry air)
    potential_temperature = 300

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 - g / (c_p * potential_temperature) * x[2]

    # pressure
    p = p0 * exner^(c_p / R)

    # temperature
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = 0
    v2 = 0

    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

@doc raw"""
    rayleigh_damping_coefficient(x, damping_rate,
                                 z_damping_start, z_damping_end,
                                 x_damping_start, x_damping_end)

Rayleigh damping coefficient of the absorbing layers used to avoid reflections at the
boundaries of a two-dimensional domain. It is the sum of the contributions of a vertical
layer, which starts at the height `z_damping_start` and ends at the top of the domain
`z_damping_end`, and of two lateral layers, which start at ``\pm`` `x_damping_start` and end
at ``\pm`` `x_damping_end`. 
"""
@inline function rayleigh_damping_coefficient(x, damping_rate,
                                              z_damping_start, z_damping_end,
                                              x_damping_start, x_damping_end)
    if x[2] <= z_damping_start
        damping_vertical = zero(damping_rate)
    else
        damping_vertical = damping_rate *
                           sinpi(0.5f0 * (x[2] - z_damping_start) /
                                 (z_damping_end - z_damping_start))^2
    end
    if x[1] < x_damping_start
        damping_right = zero(damping_rate)
    else
        damping_right = damping_rate *
                        sinpi(0.5f0 * (x[1] - x_damping_start) /
                              (x_damping_end - x_damping_start))^2
    end

    if x[1] > -x_damping_start
        damping_left = zero(damping_rate)
    else
        damping_left = damping_rate *
                       sinpi(0.5f0 * (x[1] + x_damping_start) /
                             (-x_damping_end + x_damping_start))^2
    end
    return damping_vertical + damping_right + damping_left
end

@doc raw"""
    MountainWaveSetup(; theta_0, u0, brunt_vaisala_frequency,
                      z_damping_start, z_damping_end,
                      x_damping_start, x_damping_end, damping_rate)

Setup of the mountain wave test cases, describing a uniform flow of velocity `u0` over a
stably stratified atmosphere with constant Brunt-Väisälä frequency
`brunt_vaisala_frequency` and surface potential temperature `theta_0`, i.e.
```math
\theta(z) = \theta_0 \exp(N^2 z / g), \qquad
\pi(z) = 1 + \frac{g^2}{c_p \theta_0 N^2} \big(\exp(-N^2 z / g) - 1\big).
```

The resulting object is callable with two different signatures, which are meant to be
passed to a `SemidiscretizationHyperbolic` as
- `initial_condition`, via `(x, t, equations)`, and as
- `source_terms`, via `(u, x, t, equations)`, which relaxes the solution towards the
  background state inside the absorbing layers described by `z_damping_start`,
  `z_damping_end`, `x_damping_start`, `x_damping_end`, and `damping_rate`, see
  [`rayleigh_damping_coefficient`](@ref).
"""
struct MountainWaveSetup{RealT <: Real}
    theta_0::RealT                  # background potential temperature at z = 0 in K
    u0::RealT                       # background horizontal velocity in m/s
    brunt_vaisala_frequency::RealT  # constant Brunt-Väisälä frequency in 1/s
    z_damping_start::RealT          # height where the vertical damping layer starts in m
    z_damping_end::RealT            # height where the vertical damping layer ends in m
    x_damping_start::RealT          # |x| where the lateral damping layers start in m
    x_damping_end::RealT            # |x| where the lateral damping layers end in m
    damping_rate::RealT             # Rayleigh damping rate in 1/s
end

function MountainWaveSetup(; theta_0, u0, brunt_vaisala_frequency,
                           z_damping_start, z_damping_end,
                           x_damping_start, x_damping_end, damping_rate)
    return MountainWaveSetup(promote(theta_0, u0, brunt_vaisala_frequency,
                                     z_damping_start, z_damping_end,
                                     x_damping_start, x_damping_end, damping_rate)...)
end

# Initial condition
@inline function (setup::MountainWaveSetup)(x, t,
                                            equations::AbstractCompressibleEulerEquations{2,
                                                                                          5})
    @unpack theta_0, u0, brunt_vaisala_frequency = setup
    g = equations.g
    c_p = equations.c_p
    R = equations.R    # gas constant (dry air)
    Nf = brunt_vaisala_frequency

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = 1 + g^2 / (c_p * theta_0 * Nf^2) * (exp(-Nf^2 / g * x[2]) - 1)
    # pressure
    p = equations.p_0 * exner^(c_p / R)
    potential_temperature = theta_0 * exp(Nf^2 / g * x[2])
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = u0
    v2 = 0

    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

# This signature is used for source terms, adding Rayleigh damping
# to avoid reflections at the boundaries.
@inline function (setup::MountainWaveSetup)(u, x, t,
                                            equations::CompressibleEulerPotentialTemperatureEquationsWithGravity2D)
    @unpack theta_0, u0, brunt_vaisala_frequency = setup
    @unpack z_damping_start, z_damping_end = setup
    @unpack x_damping_start, x_damping_end, damping_rate = setup

    rho, rho_v1, rho_v2, rho_theta, _ = u

    g = equations.g
    v1 = rho_v1 / rho
    theta = rho_theta / rho

    damping = -rayleigh_damping_coefficient(x, damping_rate,
                                            z_damping_start, z_damping_end,
                                            x_damping_start, x_damping_end)

    theta_b = theta_0 * exp(brunt_vaisala_frequency^2 / g * x[2])
    du2 = rho * (v1 - u0) * damping
    du3 = rho_v2 * damping
    du4 = rho * (theta - theta_b) * damping

    return SVector(zero(eltype(u)), du2, du3, du4, zero(eltype(u)))
end

@doc raw"""
    SchärMountain(; height, half_width, wavelength)

Topography of the Schär mountain test case, given by
```math
h(x) = h_c \exp\big(-(x / a)^2\big) \cos^2(\pi x / \lambda_c),
```
i.e., a train of small-scale ridges of wavelength ``\lambda_c`` = `wavelength` modulated by
a Gaussian envelope of half width ``a`` = `half_width` and peak height ``h_c`` = `height`.
The resulting object is callable as `h(x)` and is meant to be passed to
[`terrain_following_faces`](@ref).

- Christoph Schär, Daniel Leuenberger, Oliver Fuhrer, Daniel Lüthi, Claude Girard (2002)
  A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models
  Monthly Weather Review, 130(10), 2459-2480
  [DOI: 10.1175/1520-0493(2002)130](https://doi.org/10.1175/1520-0493(2002)130)
"""
struct SchärMountain{RealT <: Real}
    height::RealT      # peak height in m
    half_width::RealT  # half width of the Gaussian envelope in m
    wavelength::RealT  # wavelength of the small-scale ridges in m
end

function SchärMountain(; height, half_width, wavelength)
    return SchärMountain(promote(height, half_width, wavelength)...)
end

@inline function (mountain::SchärMountain)(x)
    @unpack height, half_width, wavelength = mountain
    return height * exp(-(x / half_width)^2) * cospi(x / wavelength)^2
end

@doc raw"""
    WitchOfAgnesi(; height, half_width)

Topography of an isolated bell-shaped mountain of peak height ``h_c`` = `height` and half
width ``a`` = `half_width`, given by the witch of Agnesi curve
```math
h(x) = \frac{h_c}{1 + (x / a)^2}.
```
The resulting object is callable as `h(x)` and is meant to be passed to
[`terrain_following_faces`](@ref). 
The flow over such a mountain is hydrostatic for
``a N / u_0 \gg 1`` and non-hydrostatic for ``a N / u_0 \sim 1``, with ``N`` the
Brunt-Väisälä frequency and ``u_0`` the background velocity, see [`MountainWaveSetup`](@ref).

- Dale R. Durran, Joseph B. Klemp (1983)
  A Compressible Model for the Simulation of Moist Mountain Waves
  Monthly Weather Review, 111(12), 2341-2361
  [DOI: 10.1175/1520-0493(1983)111](https://doi.org/10.1175/1520-0493(1983)111)
"""
struct WitchOfAgnesi{RealT <: Real}
    height::RealT      # peak height in m
    half_width::RealT  # half width in m
end

WitchOfAgnesi(; height, half_width) = WitchOfAgnesi(promote(height, half_width)...)

@inline function (mountain::WitchOfAgnesi)(x)
    @unpack height, half_width = mountain
    return height / (1 + (x / half_width)^2)
end

@doc raw"""
    terrain_following_faces(topography, x_min, x_max, z_top)

Return the tuple of the four boundary curves `(f_left, f_right, f_bottom, f_top)` of a
terrain-following domain of horizontal extent ``[x_{min}, x_{max}]`` and height `z_top`,
whose bottom boundary follows `topography`, e.g., a [`SchärMountain`](@ref) or a
[`WitchOfAgnesi`](@ref).
"""
function terrain_following_faces(topography, x_min, x_max, z_top)
    x_mid = 0.5f0 * (x_min + x_max)
    x_half_width = 0.5f0 * (x_max - x_min)
    x_of(s) = x_mid + x_half_width * s

    z_left = topography(x_min)
    z_right = topography(x_max)

    f_left(s) = SVector(x_min, z_left + (z_top - z_left) * 0.5f0 * (s + 1))
    f_right(s) = SVector(x_max, z_right + (z_top - z_right) * 0.5f0 * (s + 1))
    f_bottom(s) = SVector(x_of(s), topography(x_of(s)))
    f_top(s) = SVector(x_of(s), z_top)

    return (f_left, f_right, f_bottom, f_top)
end
